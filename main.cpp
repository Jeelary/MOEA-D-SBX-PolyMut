#include <iostream>
#include <time.h>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <limits>
#include <random>
using namespace std;
unsigned seed1 = 0;
std::mt19937 generator_norm(seed1+200);
std::normal_distribution<double> norm_dist(0.0,1.0);
double NormRand(double mu, double sigma){return norm_dist(generator_norm)*sigma + mu;}
double funk1 (double* massiv, int b)//нечётная
{
    double  sum=0, counter=0, f1=0;
        for(int j=0;j<=b;j++)
        {

            if(((j%2)==1)&&(j>=2))
            {
               counter +=1;
               sum+=(massiv[j-1]-(0.3*pow(massiv[0],2)*cos(24*3.14*massiv[0]+(4*j*3.14)/double(b))+0.6*massiv[0])*cos(6*3.14*massiv[0]+(j*3.14)/double(b)))*(massiv[j-1]-(0.3*pow(massiv[0],2)*cos(24*3.14*massiv[0]+(4*j*3.14)/double(b))+0.6*massiv[0])*cos(6*3.14*massiv[0]+(j*3.14)/double(b))); //нечетная
            }
        }
        f1=massiv[0]+(2.0/counter)*sum; //значение
       // cout<<counter<<endl;
        return f1;

}
double funk2(double* massiv,int b)//чётная
{
double  sum=0, counter=0, f2=0;
        for(int j=0;j<=b;j++)
        {

            if(((j%2)==0)&&(j>=2))
            {
               counter +=1;
               sum+=(massiv[j-1]-(0.3*pow(massiv[0],2)*cos(24*3.14*massiv[0]+(4*j*3.14)/double(b))+0.6*massiv[0])*sin(6*3.14*massiv[0]+(j*3.14)/double(b)))*(massiv[j-1]-(0.3*pow(massiv[0],2)*cos(24*3.14*massiv[0]+(4*j*3.14)/double(b))+0.6*massiv[0])*sin(6*3.14*massiv[0]+(j*3.14)/double(b))); //четная
            }
        }
        f2=1.0-sqrt(massiv[0])+(2.0/counter)*sum; //значение
        return f2;
}
   double get_betaq(const double rand, const double alpha, const double eta)
{
    double betaq = 0.0;
    if(rand <= 1.0/alpha)
        betaq = pow((rand*alpha),(1.0/(eta+1.0)));
    else
        betaq = pow((1.0/(2.0-rand*alpha)),(1.0/(eta+1.0)));
    return betaq;
}

double SBX(double eta, double y1, double y2, int j, double* Left, double* Right)
{
    if(y1 > y2)
        swap(y1,y2);
    if(y2-y1 < 10e-8)
        if(double(rand())/double(RAND_MAX) <= 0.5)
            return y1;
        else
            return y2;
    while(true)
    {
        double randval = double(rand())/double(RAND_MAX);
        if(double(rand())/double(RAND_MAX) <= 0.5)
        {
            double beta = 1.0+(2.0*(y1-Left[j])/(y2-y1));
            double alpha = 2.0-pow(beta,-(eta+1.0));
            double betaq = get_betaq(randval,alpha,eta);
            double x1_new = 0.5*((y1+y2)-betaq*(y2-y1));
            if(x1_new >= Left[j] && x1_new <= Right[j])
                return x1_new;
        }
        else
        {
            double beta = 1.0+(2.0*(Right[j]-y2)/(y2-y1));
            double alpha = 2.0-pow(beta,-(eta+1.0));
            double betaq = get_betaq(randval,alpha,eta);
            double x2_new = 0.5*((y1+y2)+betaq*(y2-y1));
            if(x2_new >= Left[j] && x2_new <= Right[j])
                return x2_new;
        }
    }
}


double PolyMut(double eta, double y, int j, double* Left, double* Right)
{

    double delta1 = (y-Left[j])/(Right[j]-Left[j]),
           delta2 = (Right[j]-y)/(Right[j]-Left[j]);
    double mut_pow = 1.0/(eta+1.0);
    double rnd = double(rand())/double(RAND_MAX), deltaq = 0.0;
    if (rnd <= 0.5)
    {
        double xy = 1.0-delta1;
        double val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta+1.0)));
        deltaq =  pow(val,mut_pow) - 1.0;
    }
    else
    {
        double xy = 1.0-delta2;
        double val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta+1.0)));
        deltaq = 1.0 - (pow(val,mut_pow));
    }
    y = y + deltaq*(Right[j]-Left[j]);
    y = std::min(Right[j], std::max(Left[j], y));
    return y;
}


int main()
{

setlocale(LC_ALL, "Russian");
int n =150, b=3, gen=200, kolFunk=2, T=7;
srand(time( NULL));
int SizeFrontPareto = 1000;
double** frontPareto = new double* [SizeFrontPareto];
double** Funk = new double *[n];//массив критериев
double** FunkY = new double *[n];//массив критериев Y
double* gteX = new double [T];
double* gteY = new double [T];
double* hi = new double [gen];
double* sredhi = new double [n];
double* memory = new double [n];
double* memory2 = new double [n];
double** mass = new double* [n];//массив случайно сгенерированных индивидов
double** massX = new double* [T];
double** massY = new double* [T];
double** EP = new double* [n];//Внешняя популяция, которая используется для хранения недоминируемых решений, найденных в ходе поиска
double** d = new double *[n]; //Евклидово расстояние между каждыми двумя точками
double** ves = new double *[n];//Равномерное распределение весовых векторов
int** B = new int *[n];
double** dmin = new double *[n];
double** Sel = new double* [n*2];
double** y = new double*[n];
double** sbx = new double*[n];
       double* L = new double[b];
       L[0] = 0;
       L[1] = -1;
       L[2] = -1;
       double* R = new double[b];
       R[0] = 1;
       R[1] = 1;
       R[2] = 1;

for (int i = 0; i < n; i++)
 {
    y[i]= new double [b];
    sbx[i]= new double [b];
    Funk[i]=new double[kolFunk];//Массив критериев
    FunkY[i]=new double[kolFunk];//Массив критериев Y
    d[i] = new double [n]; //Евклидово расстояние между каждыми двумя векторами
    dmin[i] = new double [n];//Для поиска ближайших соседей
    ves[i] = new double [kolFunk];//Весовые вектора
    B[i] = new int [n];//Массив T наиболее близких весовых векторов для каждого индивида
    mass[i]= new double[b];//Массив индивидов
    massX[i]= new double[kolFunk];
    massY[i]= new double[kolFunk];
    EP[i] = new double[b]; //Внешняя популяция (EP), которая используется для хранения недоминируемых решений, найденных в ходе поиска.
    memory2[i] = 20;
    memory[i] = 20;
 }
for (int i=0; i<gen; i++)
{
    hi[i]=0;
}
for (int i=0; i<n*2; i++)
{
    Sel[i] = new double [b];
}
for (int ig = 0; ig < 25; ig++)
{
 //ВЕСОВЫЕ ВЕКТОРА
float step1 = 0, step2 = 1;
 for(int i = 0; i < n; i++)
 {
     for(int j = 0; j < kolFunk; j++)
     {
         if(j == 0)
            ves[i][j]=step1;
         else
            ves[i][j]=step2-step1;
     }
     step1=step1+(double)1/(n-1);

 }

/*cout<<"_________________________"<<endl;
cout<<"Весовые вектора"<<endl;
cout<<"_________________________"<<endl;
 for(int i = 0; i < n; i++)
 {
     for(int j = 0; j < kolFunk; j++)
     {
         cout<<ves[i][j]<<"\t";

     }
     cout<<endl;
 }*/
//ЕВКЛИДОВЫ РАССТОЯНИЯ МЕЖДУ ДВУМЯ ВЕСОВЫМИ ВЕКТОРАМИ
double  sum=0;
for(int i = 0; i<n; i++)
 {
     int k=0;
     for(int j = 0; j<n; j++)
     {
         while(k!=kolFunk)
         {
            sum = sum + pow(ves[j][k] - ves[i][k],2);//(x2-x1)^2+(y2-y1)^2...
            k++;
         }
      d[i][j]=sqrt(sum);
      sum=0;
      k=0;
     }
 }
/*cout<<endl;
cout<<"_________________________"<<endl;
cout<<"Расстояния между векторами"<<endl;
cout<<"_________________________"<<endl;*/
for (int i = 0; i < n; i++)
{
    for (int j = 0; j < n; j++)
    {
        dmin[i][j]=d[i][j];
        //cout<<d[i][j]<<"  ";
    }
    //cout<<endl;
}
double dMin;
for (int i = 0; i < n; i++)
{

  int p=0, q=0;

while(p!=T)
{
    dMin=99999;
    for (int j = 0; j < n;j++)
    {

        if(dmin[i][j]<=dMin && dmin[i][j]!=0)
        {
            dMin=dmin[i][j];
            q=j;
        }


    }
    B[i][p]=q;
    dmin[i][q]=0;

    p++;
}

}
/*cout<<endl;
cout<<"_________________________"<<endl;
cout<<"Индексы ближайших индивидов, n="<<n<<", T="<<T<<endl;
cout<<"_________________________"<<endl;


for(int i = 0; i < n; i++)
{
    for (int j = 0; j < n; j++)
    {
       cout<<dmin[i][j]<<"  ";
    }
    cout<<endl;
}
cout<<endl;
cout<<"_________________________"<<endl;
cout<<"Индексы ближайших индивидов, n="<<n<<", T="<<T<<endl;
cout<<"_________________________"<<endl;


for(int i = 0; i < n; i++)
{
    for (int j = 0; j < T; j++)
    {
       cout<<B[i][j]<<"  ";
    }
    cout<<endl;
}*/
/*cout<<"_________________________"<<endl;
cout<<"Индивиды"<<endl;
cout<<"_________________________"<<endl;*/
  for (int i = 0; i < n; i++)
  {
        for (int j = 0; j < b; j++)
        {
         mass[i][j] = double(rand())/double(RAND_MAX);
         //cout<<mass[i][j]<<"\t";

        }
         //cout<<endl;

  }
//cout<<endl;
//cout<<"Критерии:"<<endl;
 for (int i=0;i<n;i++)
    {
        Funk[i][0]=funk1(mass[i],b);
        Funk[i][1]=funk2(mass[i],b);
        //cout<<Funk[i][0]<<" ";
        //cout<<Funk[i][1]<<endl;
    }
//cout<<endl;
//cout<<"_______________________"<<endl;
double* z = new double[kolFunk];//Референсная точка
for (int j = 0; j<kolFunk; j++)
{
    z[j]=Funk[0][j];
}
for(int j = 0; j<kolFunk; j++)
{
    for(int i = 0; i < n; i++)
    {
        if(Funk[i][j]<z[j])
          z[j]=Funk[i][j];
    }
}
/*cout<<"z= ";
for (int i = 0; i<kolFunk; i++)
{
    cout<<z[i]<<" ";
}
cout<<endl;*/
for(int g = 0; g < gen; g++)
{
int s = 0;
for(int i = 0; i < n; i++)
{
    int k = rand() % T;// Два случайных индекса из B[i]
    int l = rand() % T;
    if(k!=l)
    {
       //cout<<"_________________________ "<<endl;
       //cout<<k<<" "<<l<<endl;

       for (int j = 0; j < b; j++)
         {
            Sel[s][j] = mass[B[i][k]][j];
            Sel[s+1][j] = mass[B[i][l]][j];
         }
       s=s+2;
       int f = i * 2;
       int p = i * 2 + 1;
       hi[g] = hi[g]+memory2[i];
       memory2[i] = NormRand(memory[i],1);

              while(memory2[i]<=0)
                memory2[i] = NormRand(memory[i],1);
       for (int j = 0; j < b; j++)
         {
            sbx[i][j]=SBX(memory2[i], Sel[f][j], Sel[p][j], j, L, R);
         }


       for(int j = 0; j < b; j++)
         {
            y[i][j]=PolyMut(memory2[i],sbx[i][j],j,L,R);//новое решение y из xk и xl
         }
        FunkY[i][0]=funk1(y[i],b);//Значения критериев нового решения y
        FunkY[i][1]=funk2(y[i],b);

       for(int j = 0; j < kolFunk; j++)//Обновляем z: для каждого j=1…,m если zj > fj(y’) тогда установить zj = fj(y’)

         {
            if(z[j]>FunkY[i][j])
                z[j]=FunkY[i][j];
         }

       for(int j = 0; j < T; j++)
         {
             for(int e = 0; e < kolFunk; e++)
             {
                massX[j][e] = ves[B[i][j]][e]*abs(Funk[B[i][j]][e] - z[e]);//massX - массив значений vesj*abs(fe(xj) - ze)
                massY[j][e] = ves[B[i][j]][e]*abs(FunkY[i][e] - z[e]);//massY - массив значений vesj*abs(fe(y) - ze)
                //cout<<massX[j][e]<<"\t";
             }
             gteX[j] = massX[j][0];
             gteY[j] = massY[j][0];
             for(int e = 0; e < kolFunk; e++)
             {
                 if(gteX[j] < massX[j][e])
                    gteX[j]=massX[j][e];
                 if(gteY[j] < massY[j][e])
                    gteY[j]=massY[j][e];
             }
              if(gteY[j]<gteX[j])//если gte(y) > gte(xj), тогда установить xj=y’ и FVj=F(y’).
                {
                 int ind=0;
                 while(ind!=b)
                 {
                   mass[B[i][j]][ind]=y[i][ind];
                   ind++;
                 }
                 memory[i] = memory2[i];
                 int kf = 0;
                 while(kf!=kolFunk)
                 {
                   Funk[B[i][j]][kf]=FunkY[i][kf];
                   kf++;
                 }
                }
             }

             //cout<<endl;



       /* cout<<"X________________"<<endl;
        for(int j = 0; j < T; j++)
         {
             for(int e = 0; e < kolFunk; e++)
             {

                cout<<massX[j][e]<<"\t";
             }
             cout<<"gteX: "<<gteX[j]<<endl;

         }
         cout<<"Y________________"<<endl;
        for(int j = 0; j < T; j++)
         {
             for(int e = 0; e < kolFunk; e++)
             {

                cout<<massY[j][e]<<"\t";
             }
             cout<<"gteY: "<<gteY[j]<<endl;

         }*/




    }
    else
        i--;

}
/*cout<<"_________________________"<<endl;
cout<<"Селекция"<<endl;
cout<<"_________________________"<<endl;
for(int i=0; i<n*2; i++)
{
      for (int j = 0; j < b; j++)
            {
               cout << Sel[i][j]<<"\t";
            }
    cout << endl;
}
cout<<"_________________________"<<endl;
cout<<"Скрещивание"<<endl;
cout<<"_________________________"<<endl;
for( int i=0; i<n; i++)
{
      for (int j = 0; j < b; j++)
            {
              cout<<sbx[i][j]<<"\t";
            }
    cout << endl;
}
cout<<"________________________________"<<endl;
cout<<"Мутация"<<endl;
cout<<"________________________________"<<endl;
for(int i = 0; i < n; i++)
{
    for(int j = 0; j < b; j++)
      {
         cout<<y[i][j]<<"\t";
      }
    cout<<endl;
}
cout<<"________________________________"<<endl;
cout<<"Критерии Y"<<endl;
cout<<"________________________________"<<endl;
for(int i = 0; i < n; i++)
{
    for(int j = 0; j < kolFunk; j++)
      {
         cout<<FunkY[i][j]<<"\t";
      }
    cout<<endl;
}
cout<<"z= ";
for (int i = 0; i<kolFunk; i++)
{
    cout<<z[i]<<" ";
}
cout<<endl;
cout<<"_________________________"<<endl;
cout<<"Новые индивиды"<<endl;
cout<<"_________________________"<<endl;
  for (int i = 0; i < n; i++)
  {
        for (int j = 0; j < b; j++)
        {
         cout<<mass[i][j]<<"\t";

        }
         cout<<endl;

  }
cout<<endl;
cout<<"Новые критерии:"<<endl;
 for (int i=0;i<n;i++)
    {
        Funk[i][0]=funk1(mass[i],b);
        Funk[i][1]=funk2(mass[i],b);
        cout<<Funk[i][0]<<" ";
        cout<<Funk[i][1]<<endl;
    }
cout<<endl;
cout<<"_______________________"<<endl;*/
hi[g]=hi[g]/n;
}
double step = 0.001;
double** distP = new double* [n];
double* distPmin = new double [n];
double igd = 0, sumdist = 0;
for (int i = 0; i < SizeFrontPareto; i++)
 {
    frontPareto[i]= new double [kolFunk];
    distP[i] = new double [SizeFrontPareto];
 }

for (int i = 0; i < SizeFrontPareto; i++)
{
    for(int j = 0; j < kolFunk; j++)
    {
        frontPareto[i][j] = 0;
    }
}
for (int i = 0; i < SizeFrontPareto; i++)
{
    frontPareto[i][0] = frontPareto[i][0] + step;
    frontPareto[i][1] = 1 - sqrt(frontPareto[i][0]);
    step = step + 0.001;
}
step = 0.001;
int k=0;
sum=0;
for (int i = 0; i < n; i++)
{
    for (int j = 0; j < SizeFrontPareto; j++)
        {
           while(k!=kolFunk)
            {
               sum = sum + pow(Funk[i][k] - frontPareto[j][k], 2);//(x2-x1)^2+(y2-y1)^2...
               k++;
            }
           distP[i][j]=sqrt(sum);
           sum=0;
           k=0;
        }
}
/*for (int i = 0; i < n; i++)
{
    for (int j = 0; j < SizeFrontPareto; j++)
        {
           cout<<distP[i][j]<<"\t";
        }
        cout<<endl;
}*/
for (int i = 0; i < n; i++)
{
    dMin=99999;
    for (int j = 0; j < SizeFrontPareto; j++)
        {
          if(distP[i][j] < dMin)
          {
            distPmin[i] = distP[i][j];
            dMin = distP[i][j];
          }

        }
}
for (int i = 0; i < n; i++)
{
  sumdist = sumdist + distPmin[i];
  //cout<<distPmin[i]<<endl;
}
igd = sumdist/SizeFrontPareto;
cout<<igd<<endl;
}

ofstream index;
 index.open("h.txt");
 for(int i = 0; i < gen; i++)
{
    index<<hi[i]<<endl;;
}
index.close();

ofstream frontP;
 frontP.open("frontP.txt");
 for(int i = 0; i < SizeFrontPareto; i++)
 {
     for(int j = 0; j < kolFunk; j++)
     {
        frontP<<frontPareto[i][j]<<"\t";
     }
      frontP<<endl;
}
frontP.close();

ofstream fout;
 fout.open("mass.txt");
 for(int i=0;i<n;i++)
 {
     for(int j=0; j<b;j++)
     {
        fout<<mass[i][j]<<"\t";
     }
      fout<<endl;
}
fout.close();

ofstream funk;
funk.open("funk.txt");
for(int i=0;i<n;i++)
 {
     for(int j=0; j<2;j++)
     {
        funk<<Funk[i][j]<<"\t";
     }
     funk<<endl;
 }
 funk.close();
}












