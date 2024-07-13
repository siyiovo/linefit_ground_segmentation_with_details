# 基于Navigation框架中LiDAR地面分割算法——linefit_ground_segmentation

***该算法通过三维点云(x,y,z)降维至(d=sqrt(x^2+y^2),z)，通过判断高度z，实现激光雷达的地面分割。该算法可分为三个步骤：点云投影、直线拟合、点云分类。***

## 学习路线

### 0、前置理论要求：STL与PCL

### 1.大致读一遍论文，[详见](https://ieeexplore.ieee.org/document/5548059/figures#figures)

### 2.结合代码，看文章，[详见](https://blog.csdn.net/lovelyaiq/article/details/118826534)

***这个算法的思想比较明显，但是实现的代码比较复杂。这个内容将作为导航学习中后期的学习资料。本人学习了大概两天，主要还是代码搞得我头晕，如果对这个算法思想有一个清楚的认知，代码理解也只是时间问题。我从0到全部注释完毕，怎么也写了300多行注释，希望下一届的同学们能认真看完，这对求参数以及代码思维有比较大的帮助***

![1](https://github.com/user-attachments/assets/a28a1201-2f6b-417a-b4f6-d2e11da0e7cd)

**主要学习蓝圈里面的include与src**