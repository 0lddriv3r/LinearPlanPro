# LinearPlanPro
Linear Plan Pro
上次介绍的单纯形法大家应该发现了一个问题，在每次迭代的时候，重复进行了大量的矩阵运算，事实上我们只会用到资源向量和非基变量的系数向量。
实质上就是不断更换新基并左乘B-1，所以计算出B-1是整个算法的关键。
算法流程图如下：
