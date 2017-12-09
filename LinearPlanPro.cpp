//改进单纯形法.cpp
//-----------------------Created By 卓教授-----------------------
//说明：输入文件中问题需化为标准型(求最小值，资源向量在等式左端)
#include <iostream>
using namespace std;
#include <cmath>
#include <iomanip>
#include <fstream>
#include <ctime>

//全局变量-----------------------------------------------------------------------------------------------------------------------------------
const float M=pow(10.0,5);
int row=0,col=0,i,j,k,*R,*S,num=0,counter = 0;//counter：人工变量个数
float **A,*C,*b,*CB,**N,**tempBN,*tempb,**e,*beta,*test,*rep,*y,**B,**tempB;

//函数---------------------------------------------------------------------------------------------------------------------------------------
void newspace();   //动态开辟数组
void delspace(float **X,int n);//释放数组空间
void init(float **X,int m,int n);//初始化数组
void initI(float **X,int m);//设单位矩阵I，用以实现变换矩阵E
void output();//在屏幕输出
int maxof(float *X,int m);//求最大数
int minof(float *X,int m);//求最小正数
bool isAllNegative(float *X,int m);//判断数组是否非正数
bool isNoSolve();//判断是否无解
bool isNoLimitSolve();//判断是否为无界解
bool isMuchSolves();//判断是否有无穷多解
//float matrixMultiplication(float **B,float *b);//矩阵乘法
//------------------------------------------------------------------------------------------------------------------------------------------

int main()
{
    time_t begin,end;  //计算程序运行时间
    begin = clock();

    int subr,subc;
    float temp=0.0;

    char fileName[20],inFileName[20],outFileName[20];

    //手动输入文件名
    cout<<"Please enter the name of file(***.txt):";
    cin>>fileName;

    //定义读出文件
    strcpy(inFileName,fileName);
    strcat(inFileName,".txt");

    ifstream inFile(inFileName,ios::in);

    if (!inFile)
    {
        cout<<"Could not open the file.\n";
        exit(0);
    }

    //定义写入文件
    strcpy(outFileName,fileName);
    strcat(outFileName,"_result.txt");

    ofstream outFile(outFileName,ios::out);

    if (!outFile)
    {
        cout<<"Could not open the file.\n";
        exit(0);
    }

    //从文件中读入条件
    inFile>>row>>col;

    newspace();
    C = new float [col+row];
    b = new float [row];

    //目标函数求最小值，资源向量在等式左端。所以将b、C全部取反，即得标准化模型
    for(i=0;i<col;i++)
    {
        inFile>>C[i];
        C[i] = -C[i];
    }

    for(i=0;i<row;i++)
    {
        inFile>>b[i];
        b[i] = -b[i];
    }

    for (i=0;i<row;i++)
        for (j=0;j<col;j++)
            inFile>>A[i][j];

    //设置单位矩阵，用以实现变换矩阵E
    e = new float *[row];
    for (i=0;i<row;i++)
        e[i]=new float [row];
    initI(e,row);

    B = new float *[row];//初始可行基的逆矩阵
    for (i=0;i<row;i++)
        B[i]=new float [row];
    init(B,row,row);

    tempB = new float *[row];
    for (i=0;i<row;i++)
        tempB[i]=new float [row];
    init(tempB,row,row);

    tempb = new float [row];
    for (i=0;i<row;i++)
        tempb[i]=0;

    N = new float *[row];//非基变量系数矩阵
    for (i=0;i<col;i++)
        N[i]=new float [col];
    init(N,row,col);

    tempBN = new float *[row];//基变量的逆矩阵乘以非基变量系数矩阵
    for (i=0;i<col;i++)
        tempBN[i]=new float [col];
    init(tempBN,row,col);

    y = new float [row];//换入列的系数变换向量
    beta = new float [row];//基变量的逆矩阵乘以换入列的系数矩阵
    CB = new float [row];//基变量的价值系数
    R = new int [row];//基变量下标
    S = new int [col];//非基变量下标
    test = new float [col];//非基变量检验数
    rep = new float [row];//theta值

    for (i=0;i<row;i++)//初始化基变量下标
        R[i]=0;

    //为每一个约束条件添加人工变量，并为系数向量添上-M，并确定初始基变量下标
    for (i=0;i<row;i++)
    {
        C[col+counter] = -M;
        R[i] = col+counter;   
        counter++;
    }

    //为系数矩阵加上人工变量系数
    for(i=0;i<row;i++)
        for(j=0;j<counter;j++)
            A[i][j+col] = e[i][j];

    //确定基变量的价值系数
    for (i=0;i<row;i++)
        CB[i]=C[R[i]];

    //初始基变量全为人工变量，初始基可行解的逆矩阵为单位阵
    for (i=0;i<row;i++)
        for (j=0;j<row;j++)
            B[i][j] = e[i][j];

    //确定初始非基变量下标
    for(i=0;i<col;i++)
        S[i]=i;

    //计算非基变量系数矩阵
    for (i=0;i<row;i++)
        for (j=0;j<col;j++)
            N[i][j] = A[i][S[j]];

    //计算非基变量的检验数

    //计算基变量的逆矩阵B乘以非基变量系数矩阵N
    for (i=0;i<row;i++)
    {
        for (j=0;j<col;j++)
        {
            temp = 0;

            for (int k=0;k<row;k++)
            {
                temp += B[i][k]*N[k][j];
            }

            tempBN[i][j] = temp;
        }
    }

    for (i=0;i<col;i++)
    {
        temp = 0;

        for (j=0;j<row;j++)
        {
            temp += CB[j]*tempBN[j][i];
        }

        test[i] = C[S[i]]-temp;//test[i]对应S[i]
    }

    //计算基变量的逆矩阵B乘以资源向量b:tempb
    for (i=0;i<row;i++)
    {
        temp=0;

        for (j=0;j<row;j++)
        {
            temp += B[i][j]*b[j];
        }

        tempb[i] = temp;
    }

    cout<<"初始单纯形表:\n";
    output();

//--------------算法循环开始--------------//
    while(!isAllNegative(test,col))
    {   
        subc = maxof(test,col);   //获得检验数的最大值下标

        //计算beta向量
        cout<<"beta:";
        for (i=0;i<row;i++)
        {
            temp = 0;

            for (j=0;j<row;j++)
            {
                temp += B[i][j]*A[j][S[subc]];
            }

            beta[i] = temp;
            cout<<beta[i]<<'\t';
        }
        cout<<endl;

        //判断是否为无界解
        if (isNoLimitSolve())
        {
            outFile<<"线性方程具有无界解。\n";
            break;//若为无界解，则跳出while循环
        }//若不是无界解，则进行基变换

        //求θ值
        cout<<"rep:";
        for(i=0;i<row;i++)
        {
            if (beta[i] > 0)
            {
                rep[i] = tempb[i]/beta[i];
            }
            else
                rep[i] = -1;

            cout<<rep[i]<<'\t';
        }
        cout<<endl;

        subr=minof(rep,row);//获得θ最小值对应的下标

        num++;
        cout<<"-------------------------------------------迭代次数:"<<num<<"------------------------------------------"<<endl;

        cout<<'x'<<S[subc]+1<<"入基"<<",x"<<R[subr]+1<<"出基"<<endl;//输出出基入基情况

        //计算变换列
        temp = beta[subr];
        y[subr] = 1/temp;

        for (i=0;i<row;i++)
        {
            if (i != subr)
            {
                y[i] = -beta[i]/temp;
            }
        }

        //初始化单位阵e
        initI(e,row);

        //计算新基的变换矩阵E
        for (i=0;i<row;i++)
        {
            e[i][subr] = y[i];
        }

        //计算新基的逆矩阵
        for (i=0;i<row;i++)
        {
            for (j=0;j<row;j++)
            {
                temp = 0;

                for (k=0;k<row;k++)
                {
                    temp += e[i][k]*B[k][j];
                }

                tempB[i][j] = temp;
            }
        }

        /*for (i=0;i<row;i++)
        {
            for (j=0;j<row;j++)
            {
                for (k=0;k<row;k++)
                {
                    if (i != subr)
                        tempB[i][j] = 1*B[i][j]+y[i]*B[subr][j];
                    else
                        tempB[i][j] = y[i]*B[i][j];
                }
            }
        }*///由于变换矩阵只是单位阵换了一列，所以可以简化计算B

        cout<<"B:"<<endl;
        for (i=0;i<row;i++)
        {
            for (j=0;j<row;j++)
            {
                B[i][j] = tempB[i][j];

                cout<<B[i][j]<<'\t';
            }
            cout<<endl;
        }

        //处理下标（换基）
        temp = R[subr];
        R[subr] = S[subc];
        S[subc] = temp;

        //计算新的tempb
        for (i=0;i<row;i++)
        {
            temp = 0;

            for (j=0;j<row;j++)
            {
                temp += B[i][j]*b[j];
            }

            tempb[i] = temp;
        }

        //计算基变量的价值系数
        for (i=0;i<row;i++)
        {
            CB[i]=C[R[i]];
        }

        //计算非基变量的检验数

        //确定新的非基变量系数矩阵N
        for (i=0;i<row;i++)
        {
            for (j=0;j<col;j++)
            {
                N[i][j] = A[i][S[j]];
            }
        }

        //计算新的基变量的逆矩阵乘以非基变量的系数矩阵
        for (i=0;i<row;i++)
        {
            for (j=0;j<col;j++)
            {
                temp = 0;

                for (int k=0;k<row;k++)
                {
                    temp += B[i][k]*N[k][j];   
                }

                tempBN[i][j] = temp;
            }
        }

        cout<<"test:";
        for (i=0;i<col;i++)
        {
            temp = 0;

            for (j=0;j<row;j++)
            {
                temp += CB[j]*tempBN[j][i];
            }

            test[i] = C[S[i]]-temp;

            cout<<test[i]<<'\t';
        }
        cout<<endl;

    }//while
//-----------算法循环终止---------------//

    if (!isNoLimitSolve())
    {
        //求资源向量的值
        for (i=0;i<row;i++)
        {
            b[i] = tempb[i];
        }

        if (isNoSolve())//判断无解
        {
            outFile<<"线性方程无可行解。\n";
        }

        //存在可行解
        else
        {
            //求值
            temp=0;
            for (i=0;i<row;i++)
            {
                temp += b[i]*CB[i];
            }

            outFile<<"问题存在最优解，最优解为: \n";

            for(i=0;i<row;i++)
                outFile<<"x"<<R[i]+1<<'='<<b[i]<<'\t';

            outFile<<"\n目标函数值为：z="<<setprecision(8)<<-temp;//由于原问题是求最小值，转化为求最大值，且资源向量取反了，故目标函数值应该取反

            if (isMuchSolves())//判断是否为无穷多最优解
            {
                outFile<<"\n非基检验数中含0，问题具有无穷多最优解。\n";
            }//if
        }//else
    }//if

    cout<<"最终单纯形表:\n";
    output();

    //释放内存空间
    delspace(A,row);
    delspace(e,row);
    delspace(B,row);
    delspace(tempB,row);
    //delspace(N,row);            //bug
    //delspace(tempBN,row);        //bug
    delete[]C;
    delete[]b;
    delete[]tempb;
    delete[]beta;
    delete[]CB;
    delete[]test;
    delete[]rep;
    delete[]y;
    delete[]R;
    delete[]S;   

    cout<<"感谢使用本程序，你要求解的方程结果已存入"<<outFileName<<endl;

    end = clock();
    outFile<<"\nruntime:   "<<float(end-begin)/CLOCKS_PER_SEC<<'s'<<endl;//计算程序运行时间(文件)
    cout<<"\nruntime:   "<<float(end-begin)/CLOCKS_PER_SEC<<'s'<<endl;//计算程序运行时间(屏幕)

    return 0;
}

void newspace()   //动态开辟数组
{
    A = new float *[row];

    if (A==NULL)
    {
        cout<<"开辟数组不成功！程序结束！"<<endl;
        return ;
    }
    for (i=0;i<row;i++)
    {
        A[i]=new float[row+col];

        if(A[i]==NULL)
        {
            cout<<"开辟数组不成功！程序结束！"<<endl;
            return ;
        }
    }
}

void delspace(float **X,int m)  //释放数组空间
{
    int i;
    for(i=0;i<m;i++)
        delete []X[i];
    delete []X;
}

void init(float **X,int m,int n)
{
    for (i=0;i<m;i++)
        for (j=0;j<n;j++)
            X[i][j]=0;
}

void initI(float **X,int m)  //设置单位矩阵，用以实现变换矩阵E
{
    for (i=0;i<m;i++)
    {
        for (j=0;j<m;j++)
        {
            if (i == j)
                e[i][j] = 1;
            else
                e[i][j] = 0;
        }
    }
}

int maxof(float *X,int m)
{
    int sub = 0;
    float temp=X[0];

    for (i=1;i<m;i++)
    {
        if (X[i]>temp)
            temp = X[i];
    }

    for(i=0;i<m;i++)
        if(X[i]==temp)
            sub = i;
    return sub;
}

int minof(float *X,int m)//获得非负θ值的最小值
{
    int sub = 0;
    float temp;

    if (X[maxof(X,m)]>0)
    {
        temp = X[maxof(X,m)];
    }

       for (i=0;i<m;i++)
    {
        if (X[i]>0 && X[i]<temp)
            temp = X[i];
    }

    for(i=0; i<m; i++)
        if(X[i] == temp)
            sub = i;

    return sub;
}

/*void matrixMultiplication(float **B,float *b,float **result)
{
    for (i=0;i<row;i++)
    {
        for (j=0;j<row;j++)
        {
            temp = 0;
            for (int k=0;k<row;k++)
            {
                temp += B[i][k]*b[k][j];
                result[i][j] = temp;
            }
        }
    }
}*/

bool isAllNegative(float *X,int m)
{
    int flag = 0;
    for (i=0;i<m;i++)
        if(X[i] > 0)
            flag += 1;
    if (flag == 0)
        return true;
    else
        return false;
}

//判断是否无解
bool isNoSolve()
{
    int flag = 0;
    for(i=0;i<row;i++)
    {
        if(R[i]>col && b[i]!=0)//基变量中是否含有非零的人工变量
            flag += 1;
    }
    if (flag == 0)
        return false;
    else
        return true;
}

bool isNoLimitSolve()//判断是否为无界解
{
    if (isAllNegative(test,col))
        return 0;
    else
    {
        int flag=0,bol=0;

        for (j=0;j<row;j++)
        {
            if (beta[j]>0)//换入变量在单纯性表中对应的列如果全部<=0，则为无界解
                flag += 1;
        }

        if (flag != 0)
            bol += 1;

        if (bol == 0)
            return 1;
        else
            return 0;
    }
}

bool isMuchSolves()//判断是否有无穷多解（非基变量的检验数中如果存在一个等于0，则有无穷多最优解）
{
    for (i=0;i<col;i++)
    {
        if (fabs(test[i]) < 1.0e-7)//考虑计算机的舍入误差，即：fabs(test[i]) < 1.0e-7
            return 1;
    }
    return 0;
}

void output()
{
    cout<<"-------------------------------------------------------------------------------\n";
    cout<<"                ";
    for (i=0;i<col+counter;i++)
        cout<<C[i]<<'\t';
    cout<<endl;
    for(i=0;i<row;i++)
    {
        cout<<CB[i]<<"   "<<b[i]<<"   "<<'\t';
        for(j=0;j<col+counter;j++)
            cout<<A[i][j]<<'\t';
        cout<<endl;
    }
    cout<<"                ";
    for(i=0;i<col;i++)
        cout<<test[i]<<'\t';
    cout<<"\n-------------------------------------------------------------------------------\n";

}
