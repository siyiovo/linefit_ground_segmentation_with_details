#include "ground_segmentation/segment.h"

Segment::Segment(const unsigned int& n_bins,
                 const double& min_slope,
                 const double& max_slope,
                 const double& max_error,
                 const double& long_threshold,
                 const double& max_long_height,
                 const double& max_start_height,
                 const double& sensor_height) :
                 bins_(n_bins),
                 min_slope_(min_slope),
                 max_slope_(max_slope),
                 max_error_(max_error),
                 long_threshold_(long_threshold),
                 max_long_height_(max_long_height),
                 max_start_height_(max_start_height),
                 sensor_height_(sensor_height){}

void Segment::fitSegmentLines() {
  // Find first point.首先找到该segment下第一个含有(d,z)的bin，详见insertionThread
  auto line_start = bins_.begin();  
  while (!line_start->hasPoint()) {  //遍历bin下的所有点，直到找到第一个含有2D坐标(d,z)的bin
    ++line_start;
    // Stop if we reached last point.
    if (line_start == bins_.end()) return;  //若这个segment下bin无点(d,z)，返回
  }
  // Fill lines.
  bool is_long_line = false;  //该直线长度是否超出long_threshold_这个阈值
  double cur_ground_height = -sensor_height_;  //以LiDAR所在的x0y平面为基准面，底盘面的高度就是负值
  std::list<Bin::MinZPoint> current_line_points(1, line_start->getMinZPoint());  //储存刚刚找到的bin，这是要拟合直线的起点
  LocalLine cur_line = std::make_pair(0,0);  //要拟合直线的斜率k与偏移b
  for (auto line_iter = line_start+1; line_iter != bins_.end(); ++line_iter) {  //line_iter始于bin_2(bin_1 + 1)
    if (line_iter->hasPoint()) {  //如果该bin有(d,z)
      Bin::MinZPoint cur_point = line_iter->getMinZPoint();  //获取该bin的(d,z)
      if (cur_point.d - current_line_points.back().d > long_threshold_) is_long_line = true;  //如果该bin.d - 起点的模长d是大于segments_设置的长直线阈值，判断为长直线
      if (current_line_points.size() >= 2) {  //注：第一次迭代无法进入这个判断，因为此时current_line_points内只含有bin_1，bin_2还未判断是否与bin_1在同一水平面
        // Get expected z value to possibly reject far away points.假设已经判断已经拟合完第一条直线，且第一条直线与第二条直线斜率不同
        double expected_z = std::numeric_limits<double>::max();  //double类型的极大值，为1.79769e+308，保护expected_z
        if (is_long_line && current_line_points.size() > 2) {  //如果第三个点.d-第二个点.d是长直线
          expected_z = cur_line.first * cur_point.d + cur_line.second;  //在拟合到一条真正的直线之前，expected_z即为0. z = k * d + b
        }
        current_line_points.push_back(cur_point);  //添加第三个点，现在当前拟合直线的点集里面有3个点了
        cur_line = fitLocalLine(current_line_points);  //直线拟合，求斜率k与偏移b
        const double error = getMaxError(current_line_points, cur_line);  //求z的理论值与最大值的最大偏差
        // Check if not a good line.先不判断第三个点是否满足拟合，而是先求出k、b等值，再拿这些值来判断这个点是否舍去，非常好想法，爱来自中国
        if (error > max_error_ ||
            std::fabs(cur_line.first) > max_slope_ ||
            (current_line_points.size() > 2 && std::fabs(cur_line.first) < min_slope_) ||
            is_long_line && std::fabs(expected_z - cur_point.z) > max_long_height_) {  //如果求出来的值大于/小于相应阈值的话，证明第三个点不符合现在拟合的这条曲线
          // Add line until previous point as ground.  
          current_line_points.pop_back();  //删掉第三个点
          // Don't let lines with 2 base points through.
          if (current_line_points.size() >= 3) {  //如果第三个点保留，进入判断
            const LocalLine new_line = fitLocalLine(current_line_points);  //再拟合，第三个点还保留的跟之前一致
            lines_.push_back(localLineToLine(new_line, current_line_points));  //把线以起点与终点的形式记下来，(d,z[i])
            cur_ground_height = new_line.first * current_line_points.back().d + new_line.second;  //重新定义底盘面高度
          }
          // Start new line.
          is_long_line = false;
          current_line_points.erase(current_line_points.begin(), --current_line_points.end());  //删掉前两个点
          --line_iter;  //下一次迭代从最后一个（该示例里的第三个）点开始，为了连接相邻的两条直线
        }
        // Good line, continue.
        else { }
      }
      else {
        // Not enough points.一般来说这是一条直线中第一次迭代要判断的
        if (cur_point.d - current_line_points.back().d < long_threshold_ &&
            std::fabs(current_line_points.back().z - cur_ground_height) < max_start_height_) {  //如果这玩意不是长直线，并且以LiDAR为基准面时，第二个点的高度小于拟合直线的最大高度（第一个点到第二个点的高度差小于某个常数）
          // Add point if valid.添加该直线的第二个点，这个时候可以进入上一个判断中
          current_line_points.push_back(cur_point);
        }
        else {  //如果这直线已经足够长了，或者说第一个点与第二个点高度差过大
          // Start new line.重新找起点
          current_line_points.clear();
          current_line_points.push_back(cur_point);  //舍弃掉bin_1，这样子下一次迭代的时候就会将第二个点设为起点
        }
      }
    }
  }
  // Add last line.很明显，当处理完最后一条直线后，没有进行该直线内点的拟合
  if (current_line_points.size() > 2) {
    const LocalLine new_line = fitLocalLine(current_line_points);
    lines_.push_back(localLineToLine(new_line, current_line_points));  //至此，同一个segment内的所有bin内所有能拟合的直线全部拟合完成，接下来开始点云分类，详见assignCluster
  }
}

Segment::Line Segment::localLineToLine(const LocalLine& local_line,
                                       const std::list<Bin::MinZPoint>& line_points) {
  Line line;
  const double first_d = line_points.front().d;
  const double second_d = line_points.back().d;
  const double first_z = local_line.first * first_d + local_line.second;
  const double second_z = local_line.first * second_d + local_line.second;
  line.first.z = first_z;  //起点.拟合后的z
  line.first.d = first_d;  //起点.d
  line.second.z = second_z;  //终点.拟合后的z
  line.second.d = second_d;  //终点.d
  return line;
}

double Segment::verticalDistanceToLine(const double &d, const double &z) {
  static const double kMargin = 0.1;
  double distance = -1;
  for (auto it = lines_.begin(); it != lines_.end(); ++it) {  //计算segment[i]的所有直线的投影误差
    if (it->first.d - kMargin < d && it->second.d + kMargin > d) {  /*point_2d(z,d)代表该线程里头算有点云的2D坐标
                                                                    kMargin余量用于保证当前迭代的点云2D坐标在所有直线中的任意一条里头
                                                                    所以kMargin余量的值比较重要，不然该点可能不属于任何一条直线*/
      const double delta_z = it->second.z - it->first.z;  //终点.z-起点.z
      const double delta_d = it->second.d - it->first.d;  //终点.d-起点.d
      const double expected_z = (d - it->first.d)/delta_d *delta_z + it->first.z;  //已知斜率k=delta_z/delta_d，易得直线的截距式为expected_z = k ( d - 起点.d ) + 起点.z，与这行代码表达式一致
      distance = std::fabs(z - expected_z);  //计算投影误差
    }
  }
  return distance;
}

double Segment::getMeanError(const std::list<Bin::MinZPoint> &points, const LocalLine &line) {
  double error_sum = 0;
  for (auto it = points.begin(); it != points.end(); ++it) {
    const double residual = (line.first * it->d + line.second) - it->z;
    error_sum += residual * residual;
  }
  return error_sum/points.size();
}

double Segment::getMaxError(const std::list<Bin::MinZPoint> &points, const LocalLine &line) {
  double max_error = 0;
  for (auto it = points.begin(); it != points.end(); ++it) {
    const double residual = (line.first * it->d + line.second) - it->z;  //residual残差，delta_z=z[i]-z，拟合后点的高度z[i] = k * d[i] + b
    const double error = residual * residual;  //偏差
    if (error > max_error) max_error = error;
  }
  return max_error;
}

Segment::LocalLine Segment::fitLocalLine(const std::list<Bin::MinZPoint> &points) {
  const unsigned int n_points = points.size();  //假设现在在第二次迭代里，已知该直线有3个点
  Eigen::MatrixXd X(n_points, 2);  //矩阵X是3x2维矩阵
  Eigen::VectorXd Y(n_points);  //矩阵Y是3维向量
  unsigned int counter = 0;
  for (auto iter = points.begin(); iter != points.end(); ++iter) {  //传参
    X(counter, 0) = iter->d;  //把各点的模长d赋值到X第一列里
    X(counter, 1) = 1;  //X第二列就赋值1
    Y(counter) = iter->z;  //Y依次赋值各地的高度z（LiDAR所在平面是基准面）
    ++counter;
  }
  Eigen::VectorXd result = X.colPivHouseholderQr().solve(Y);  //相当于X*A=Y，求A
  LocalLine line_result;  //返回拟合完的k与b
  line_result.first = result(0);
  line_result.second = result(1);
  return line_result;
}

bool Segment::getLines(std::list<Line> *lines) {
  if (lines_.empty()) {
    return false;
  }
  else {
    *lines = lines_;
    return true;
  }
}
