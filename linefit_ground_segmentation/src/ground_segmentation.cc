#include "ground_segmentation/ground_segmentation.h"

#include <chrono>
#include <cmath>
#include <list>
#include <memory>
#include <thread>

using namespace std::chrono_literals;


GroundSegmentation::GroundSegmentation(const GroundSegmentationParams& params) :
    params_(params),
    segments_(params.n_segments, Segment(params.n_bins,
                                         params.min_slope,
                                         params.max_slope,
                                         params.max_error_square,
                                         params.long_threshold,
                                         params.max_long_height,
                                         params.max_start_height,
                                         params.sensor_height)) {
  if (params.visualize) viewer_.reset(new Viewer());
}

void GroundSegmentation::segment(const PointCloud& cloud, std::vector<int>* segmentation) {
  std::cout << "Segmenting cloud with " << cloud.size() << " points...\n";
  std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
  segmentation->clear();  //清除容器内的点云元素，而不释放点云内存
  segmentation->resize(cloud.size(), 0);  //初始化容器大小，通过PCL1给的点云数量resize
  bin_index_.resize(cloud.size());  //初始化(segment[i],bin[i])索引的数量
  segment_coordinates_.resize(cloud.size());  //初始化每个点云在其对应的segment的2D坐标(d,z)的数量
  resetSegments();  //初始化参数segments_

  insertPoints(cloud);  //点云投影
  std::list<PointLine> lines;
  if (params_.visualize) {  //如果点云可视化，就获取拟合线的数量，这其中包含直线拟合
    getLines(&lines);  //直线拟合
  }
  else {
    getLines(NULL);
  }
  assignCluster(segmentation);  //点云分类

//到这里，地面分割算法就大概结束了，以下是拟合完直线的可视化，大家自行看看。总的来说，这个算法不难，但是代码比较复杂，如果看得很崩溃（我当时没有这些注释），就先玩玩米家游戏放松一下~

  size_t n_ground = 0;
  for (const auto seg: *segmentation) {
    n_ground += seg;
  }

  if (params_.visualize) {
    // Visualize.
    PointCloud::Ptr obstacle_cloud = boost::make_shared<PointCloud>();
    obstacle_cloud->reserve(segmentation->size() - n_ground);
    // Get cloud of ground points.
    PointCloud::Ptr ground_cloud = boost::make_shared<PointCloud>();
    ground_cloud->reserve(n_ground);
    for (size_t i = 0; i < cloud.size(); ++i) {
      if (segmentation->at(i) == 1) ground_cloud->push_back(cloud[i]);
      else obstacle_cloud->push_back(cloud[i]);
    }
    PointCloud::Ptr min_cloud = boost::make_shared<PointCloud>();
    getMinZPointCloud(min_cloud.get());
    viewer_->visualize(lines, min_cloud, ground_cloud, obstacle_cloud);
  }
  std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> fp_ms = end - start;
  std::cout << "Done! Took " << fp_ms.count() << "ms\n";
}

void GroundSegmentation::getLines(std::list<PointLine> *lines) {
  std::mutex line_mutex;
  std::vector<std::thread> thread_vec(params_.n_threads);  //依然是分n_threads个线程处理
  unsigned int i;
  for (i = 0; i < params_.n_threads; ++i) {
    const unsigned int start_index = params_.n_segments / params_.n_threads * i;
    const unsigned int end_index = params_.n_segments / params_.n_threads * (i+1);
    thread_vec[i] = std::thread(&GroundSegmentation::lineFitThread, this,
                                start_index, end_index, lines, &line_mutex);  //储存这个线程所有拟合好的直线，详见lineFitThread
  }
  for (auto it = thread_vec.begin(); it != thread_vec.end(); ++it) {
    it->join();
  }
}

void GroundSegmentation::lineFitThread(const unsigned int start_index,
                                       const unsigned int end_index,
                                       std::list<PointLine> *lines, std::mutex* lines_mutex) {
  const bool visualize = lines;
  const double seg_step = 2*M_PI / params_.n_segments;
  double angle = -M_PI + seg_step/2 + seg_step * start_index;
  for (unsigned int i = start_index; i < end_index; ++i) {
    segments_[i].fitSegmentLines();  //直线拟合
    // Convert lines to 3d if we want to.如果想直线可视化到RVIZ，就将起点、终点投影回3D坐标系中
    if (visualize) {
      std::list<Segment::Line> segment_lines;
      segments_[i].getLines(&segment_lines);
      for (auto line_iter = segment_lines.begin(); line_iter != segment_lines.end(); ++line_iter) {
        const pcl::PointXYZ start = minZPointTo3d(line_iter->first, angle);  //将二维点投影回三维点
        const pcl::PointXYZ end = minZPointTo3d(line_iter->second, angle);
        lines_mutex->lock();
        lines->emplace_back(start, end);
        lines_mutex->unlock();
      }

      angle += seg_step;
    }
  }
}

void GroundSegmentation::getMinZPointCloud(PointCloud* cloud) {
  cloud->reserve(params_.n_segments * params_.n_bins);
  const double seg_step = 2*M_PI / params_.n_segments;
  double angle = -M_PI + seg_step/2;
  for (auto seg_iter = segments_.begin(); seg_iter != segments_.end(); ++seg_iter) {
    for (auto bin_iter = seg_iter->begin(); bin_iter != seg_iter->end(); ++bin_iter) {
      const pcl::PointXYZ min = minZPointTo3d(bin_iter->getMinZPoint(), angle);
      cloud->push_back(min);
    }

    angle += seg_step;
  }
}

void GroundSegmentation::resetSegments() {
  segments_ = std::vector<Segment>(params_.n_segments, Segment(params_.n_bins,
                                                               params_.min_slope,
                                                               params_.max_slope,
                                                               params_.max_error_square,
                                                               params_.long_threshold,
                                                               params_.max_long_height,
                                                               params_.max_start_height,
                                                               params_.sensor_height));
}

pcl::PointXYZ GroundSegmentation::minZPointTo3d(const Bin::MinZPoint &min_z_point,
                                                const double &angle) {
  pcl::PointXYZ point;
  point.x = cos(angle) * min_z_point.d;
  point.y = sin(angle) * min_z_point.d;
  point.z = min_z_point.z;
  return point;
}

void GroundSegmentation::assignCluster(std::vector<int>* segmentation) {
  std::vector<std::thread> thread_vec(params_.n_threads);  //将线程储存为一个n_threads维向量
  const size_t cloud_size = segmentation->size();  //获取这一帧的点云数量
  for (unsigned int i = 0; i < params_.n_threads; ++i) {  //依然是按n_threads个线程依次处理

  //下面操作与各线程的点云投影一样
    const unsigned int start_index = cloud_size / params_.n_threads * i;  
    const unsigned int end_index = cloud_size / params_.n_threads * (i+1);
    thread_vec[i] = std::thread(&GroundSegmentation::assignClusterThread, this,
                                start_index, end_index, segmentation);
  }
  for (auto it = thread_vec.begin(); it != thread_vec.end(); ++it) {
    it->join();  //等待当前线程处理完毕
  }
}

void GroundSegmentation::assignClusterThread(const unsigned int &start_index,
                                             const unsigned int &end_index,
                                             std::vector<int> *segmentation) {
  const double segment_step = 2*M_PI/params_.n_segments;
  for (unsigned int i = start_index; i < end_index; ++i) {
    Bin::MinZPoint point_2d = segment_coordinates_[i];  //储存该线程下每个点云的2D坐标，详见insertionThread
    const int segment_index = bin_index_[i].first;  //segment索引
    if (segment_index >= 0) {
      double dist = segments_[segment_index].verticalDistanceToLine(point_2d.d, point_2d.z);  //求在每个segment下，每个点云的投影误差
      // Search neighboring segments.
      int steps = 1;
      while (dist < 0 && steps * segment_step < params_.line_search_angle) {  //当这个点真的不在他所在的semgment所有拟合曲线内时，分别求这个点云在前后两个segment的投影误差
        // Fix indices that are out of bounds.
        int index_1 = segment_index + steps;  //后一个segment
        while (index_1 >= params_.n_segments) index_1 -= params_.n_segments;
        int index_2 = segment_index - steps;  //前一个segment
        while (index_2 < 0) index_2 += params_.n_segments;
        // Get distance to neighboring lines.
        const double dist_1 = segments_[index_1].verticalDistanceToLine(point_2d.d, point_2d.z);
        const double dist_2 = segments_[index_2].verticalDistanceToLine(point_2d.d, point_2d.z);
        if (dist_1 >= 0) {
          dist = dist_1;
        }
        if (dist_2 >= 0) {
          // Select smaller distance if both segments return a valid distance.
          if (dist < 0 || dist_2 < dist) {
            dist = dist_2;
          }
        }
        ++steps;
      }
      if (dist < params_.max_dist_to_line && dist != -1) {  //如果这个点云的投影误差小于阈值，那就把它归为地面点，否则该点云归为非地面点
        segmentation->at(i) = 1;  //归为地面点
      }
    }
  }
}

void GroundSegmentation::getMinZPoints(PointCloud* out_cloud) {
  const double seg_step = 2*M_PI / params_.n_segments;
  const double bin_step = (sqrt(params_.r_max_square) - sqrt(params_.r_min_square))
      / params_.n_bins;
  const double r_min = sqrt(params_.r_min_square);
  double angle = -M_PI + seg_step/2;
  for (auto seg_iter = segments_.begin(); seg_iter != segments_.end(); ++seg_iter) {
    double dist = r_min + bin_step/2;
    for (auto bin_iter = seg_iter->begin(); bin_iter != seg_iter->end(); ++bin_iter) {
      pcl::PointXYZ point;
      if (bin_iter->hasPoint()) {
        Bin::MinZPoint min_z_point(bin_iter->getMinZPoint());
        point.x = cos(angle) * min_z_point.d;
        point.y = sin(angle) * min_z_point.d;
        point.z = min_z_point.z;

        out_cloud->push_back(point);
      }
      dist += bin_step;
    }
    angle += seg_step;
  }
}

void GroundSegmentation::insertPoints(const PointCloud& cloud) {
  std::vector<std::thread> threads(params_.n_threads);  //首先将所有点云分成n_threads份
  const size_t points_per_thread = cloud.size() / params_.n_threads;  //一个线程要处理的点云数量
  // Launch threads.
  for (unsigned int i = 0; i < params_.n_threads - 1; ++i) {
    const size_t start_index = i * points_per_thread;  //eg:[0~10000]是第一个线程要处理的点云，其中0是start_index，10000是end_index，第二个线程[10001~20000]同理
    const size_t end_index = (i+1) * points_per_thread;  //注：原作者第一次的提交是points_per_thread - 1，这样会漏掉每个线程的最后一个点云
    threads[i] = std::thread(&GroundSegmentation::insertionThread, this,
                             cloud, start_index, end_index);  //分好每个线程的任务后，开始点云投影
  }
  // Launch last thread which might have more points than others.
  const size_t start_index = (params_.n_threads - 1) * points_per_thread;  //处理最后一个线程的点云，要单独处理是因为这个线程要处理点云的数量与别的线程不一样
  const size_t end_index = cloud.size();  //同上注释
  threads[params_.n_threads - 1] =
      std::thread(&GroundSegmentation::insertionThread, this, cloud, start_index, end_index);  //同上
  // Wait for threads to finish.
  for (auto it = threads.begin(); it != threads.end(); ++it) {  //一个一个线程来
    it->join();  //直到该线程处理完，该线程环境都是block状态，防止串数据
  }
}

void GroundSegmentation::insertionThread(const PointCloud& cloud,
                                         const size_t start_index,
                                         const size_t end_index) {
  const double segment_step = 2*M_PI / params_.n_segments;  //每个segment的角度，同时也是步长
  const double bin_step = (sqrt(params_.r_max_square) - sqrt(params_.r_min_square))  //每个segment下bin[i]的正步长
      / params_.n_bins;
  const double r_min = sqrt(params_.r_min_square);
  for (unsigned int i = start_index; i < end_index; ++i) {
    pcl::PointXYZ point(cloud[i]);  //储存每个点云的3D坐标(x,y,z)
    const double range_square = point.x * point.x + point.y * point.y;  
    const double range = sqrt(range_square);  //每个点云的模长
    if (range_square < params_.r_max_square && range_square > params_.r_min_square) {  //在范围内的点云
      const double angle = std::atan2(point.y, point.x);  //每个点云的角度
      const unsigned int bin_index = (range - r_min) / bin_step;  //每个segment下每个bin的索引，即求该点云在哪个bin
      const unsigned int segment_index = (angle + M_PI) / segment_step;  //segment[i]的索引，即求该点云在哪个segment。由于atan2求出来的值域为(-pi,pi)，要把点云的角度范围设在(0,2pi)，即+ pi
      //注意segment_index_clamped，bug fix原因是因为第一次提交中没有考虑到一个点云的角度恰好为pi，此时该类点云在segment[360]处，这类点云应该算到第一个segment处，即segment_index为0
      const unsigned int segment_index_clamped = segment_index == params_.n_segments ? 0 : segment_index;
      segments_[segment_index_clamped][bin_index].addPoint(range, point.z);  //储存2D最小坐标值(d,z)到相应的索引标签里，注：这个不用来储存点云2D坐标
      bin_index_[i] = std::make_pair(segment_index_clamped, bin_index);  //给每个点云储存好索引标签，即该点云在哪个(segment,bin)
    }
    else {  //不在范围的点云
      bin_index_[i] = std::make_pair<int, int>(-1, -1);  //不在ROI里的点云就不管它，default:-1
    }
    segment_coordinates_[i] = Bin::MinZPoint(range, point.z);  /*储存每个点云的2D坐标
                                                               这个函数处理了该帧下3D->2D的点云投影，同时也划分好了区域
                                                               注：没有给点云在线程层面添加索引，线程索引只是标记线程
                                                               只有bin_index_是每个点云在当前线程的索引*/
  }
}
