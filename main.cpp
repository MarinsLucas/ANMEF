/* TO COMPILE:
 g++ *.cpp *.h -o ./output.exe
 
 */
#include<vector>
#include<cmath>
#include<iostream>
#include<eigen3/Eigen/Dense>
#include<fstream>

#define f long double

#define DIRICHLET 0
#define NEUMANN 1

#define EPSLON 0.0001


using namespace std;
using namespace Eigen;

struct contourCondition
{
    f value;
    int type;
};

contourCondition create_contourCondition(f value, int type)
{
    contourCondition cc;
    cc.value = value;
    cc.type = type;
    return cc;
}

f cospix(f x)
{
    return cos(M_PI*x);
}

f constante_1(f x)
{
    return 1; 
}

f pi2sinPiX(f x)
{
    return M_PI*M_PI*sin(M_PI*x);
}

f sinpix(f x)
{
    return sin(M_PI * x);
}

f picospix(f x)
{
    return M_PI*cos(M_PI*x);
}

f mpicospix(f x)
{
    return - M_PI*cos(M_PI*x);
}

f pi2cospix(f x)
{
    return M_PI*M_PI*cos(M_PI*x);
}
f solexata2q(f x)
{
    f epslon = EPSLON; 
    f c2 = (pow(M_E, -1.0/sqrt(epslon)) - 1.0)/(pow(M_E, 1.0/sqrt(epslon)) - pow(M_E, -1.0/sqrt(epslon)));
    f c1 = - 1.0 - c2;
    return c1*pow(M_E, -x/sqrt(epslon)) + c2*pow(M_E, x/sqrt(epslon)) + 1.0;
}

void printvector(vector<f> v)
{
    for(int i = 0; i<v.size(); i++)
    {
        cout<<v[i]<<",";
    }
    cout<<endl;
}

void printmatrix(MatrixXd m)
{
    ofstream file("matrix.csv");
    if(!file)
    {
        cerr<<"Erro ao abrir arquivo de saída"<<endl;
        return;
    }
    file<<m<<endl;
}


vector<f> gauss_weights(int nint)
{
    vector<f> w(nint, 0);
    switch (nint)
    {
    case 2:
        w[0] = 1.0;
        w[1] = 1.0;
        break;

    case 3:
        w[0] = 5/9.0;
        w[1] = 8/9.0;
        w[2] = 5/9.0;
        break;
    
    case 4:
        w[0] = (18.0 - sqrt(30.0))/36.0;
        w[1] = (18.0 + sqrt(30.0))/36.0;
        w[2] = (18.0 + sqrt(30.0))/36.0;
        w[3] = (18.0 - sqrt(30.0))/36.0;
        break;
    
    case 5:
        w[0] = (322 - 13*sqrt(70))/900.0;
        w[1] = (322 + 13*sqrt(70))/900.0;
        w[2] = 128.0/225.0;
        w[3] = (322 + 13*sqrt(70))/900.0;
        w[4] = (322 - 13*sqrt(70))/900.0;
        break;

    default:
        cout<<"gauss_weights error: não implementei esse número de pontos de integração ainda"<<endl;
        break;
    }
    return w;
}

f*** create_shg(int nen, int nint)
{
    f*** shg = new f**[2];
    shg[0] = new f*[nen];
    shg[1] =  new f*[nen];

    vector<f> pt;
    switch (nint)
    {
    case 2:
        /* pt.push_back(-1/sqrt(3.0));
        pt.push_back(1/sqrt(3.0)); */
        pt.push_back(-1.0);
        pt.push_back(1.0);
        break;
    
    case 3:
        pt.push_back(-sqrt(3.0/5.0));
        pt.push_back(0.0);
        pt.push_back(sqrt(3.0/5.0));
        break;
    
    case 4:
        pt.push_back(-sqrt((3.0/7.0) + (2.0/7.0)*sqrt(6.0/5.0)));
        pt.push_back(-sqrt((3.0/7) - (2.0/7.0)*sqrt(6.0/5.0)));
        pt.push_back(sqrt((3.0/7) - (2.0/7.0)*sqrt(6.0/5.0)));
        pt.push_back(sqrt((3.0/7) + (2.0/7.0)*sqrt(6.0/5.0)));
        break;

    case 5:
        pt.push_back((-1.0/3.0)*sqrt(5.0 +  2.0*sqrt(10.0/7.0)));
        pt.push_back((-1.0/3.0)*sqrt(5.0 -  2.0*sqrt(10.0/7.0)));
        pt.push_back(0.0);
        pt.push_back((1.0/3.0)*sqrt(5.0 -  2.0*sqrt(10.0/7.0)));
        pt.push_back((1.0/3.0)*sqrt(5.0 +  2.0*sqrt(10.0/7.0)));
        break;

    default:
        cout<<"shg error: não implementei esse número de pontos de integração ainda"<<endl;

        break;
    }

    
    for(int j=0;j<nen; j++)
    {
        shg[0][j] = new f[nint];
        shg[1][j] = new f[nint];
    }
    
    //Preenchimento do shg
    for(int l = 0; l<nint; l++)
    {
        f t = pt[l];
        switch (nen)
        {
        case 2: 
            shg[0][0][l] = (1.0 - t)/2.0;
            shg[0][1][l] = (1.0 + t)/2.0;
            shg[1][0][l] = -1.0/2.0;
            shg[1][1][l] = 1.0/2.0;
            break;
        
        case 3:               
            shg[0][0][l] = (1/2.0)*t*(t-1.0);
            shg[0][1][l] = -((t-1.0)*(t+1.0));
            shg[0][2][l] = (1/2.0)*t*(t+1);
            shg[1][0][l] = (2.0*t - 1)/2.0;
            shg[1][1][l] = -(2.0*t);
            shg[1][2][l] = (2.0*t + 1)/2.0;
            break;

        case 4:
            shg[0][0][l] = -(9.0/16.0)*(t + 1.0/3.0)*(t - 1.0/3.0)*(t-1.0);
            shg[0][1][l] = (27.0/16.0)*(t + 1)*(t - 1.0/3.0)* (t-1);
            shg[0][2][l] =  -(27.0/16.0)*(t + 1.0)*(t+ 1.0/3.0)*(t - 1.0);
            shg[0][3][l] = (9.0/16.0)*(t + 1.0)*(t+1.0/3.0)*(t - 1.0/3.0);

            shg[1][0][l] = -(9.0/16.0)*(3.0*pow(t,2) - 2.0*t - 1.0/9.0);
            shg[1][1][l] = (27.0/16.0)*(3.0*pow(t, 2) - 2.0*t/3.0 - 1);
            shg[1][2][l] =  (-27.0/16.0)*(3.0*pow(t, 2) + 2.0*t/3.0 -1);
            shg[1][3][l] =  (9.0/16.0) * (3.0*pow(t, 2) + 2.0*t - 1.0/9.0);
            break;

        case 5:
            shg[0][0][l] = (t + (1.0/2.0))*t*(t - (1.0/2.0))*(t - 1.0)*(2.0/3.0);
            shg[0][1][l] = (t + 1.0)*t*(t - (1.0/2.0))*(t - 1.0)*(-8.0/3.0);
            shg[0][2][l] =  (t + 1.0)*(t + (1.0/2.0))*(t - (1.0/2.0))*(t - 1.0)*4.0;
            shg[0][3][l] = (t + 1.0)*(t + (1.0/2.0))*t*(t - 1.0)*(-8.0/3.0);
            shg[0][4][l] = (t + 1.0)*(t + (1.0/2.0))*t*(t - (1.0/2.0))*(2.0/3.0);

            shg[1][0][l] = ((16.0*pow(t, 3.0)) - (12.0*pow(t, 2.0)) - (2.0*t) + 1.0)/6.0;
            shg[1][1][l] = -4.0*((8.0*pow(t,3.0)) - (3.0*pow(t,2.0)) - (4.0*t) + 1.0)/3.0;
            shg[1][2][l] =  (16.0*pow(t, 3)) - (10.0*t);
            shg[1][3][l] = -4*((8.0*pow(t,3.0)) + (3.0*pow(t,2.0)) - (4.0*t) - 1.0)/3.0;
            shg[1][4][l] = ((16.0*pow(t,3.0)) + (12.0*pow(t,2.0)) - (2*t) - 1.0)/6.0;
            break;

        default:
            cout<<"Não definida função shg para esse número de pontos de integração"<<endl;
            break;
        }
    }

    return shg; 
}

f errul2(int nel, int nint, int nen, VectorXd u, f (*solexata)(f), f h, vector<f> x, short t)
{
    //Inicialize weight vector
    vector<f> w = gauss_weights(nint);

    //Inicializando shg
    f*** shg = create_shg(nen, nint);

    f erl2 = 0.0;
    int index = 0;
    for(int n =0; n<nel; n++)
    {
        f eru = 0;
        for(int l = 0; l < nint; l++)
        {
            f uh = 0;
            f xx = 0;

            for(int i=0; i<nen; i++)
            {
                xx = xx + shg[t][i][l] *(x[index + i]);
                if(t==0)
                {
                    uh = uh + shg[t][i][l] * u(index + i); //Tenho que selecionar o u da mesma forma do de baixo
                }else{
                    uh = uh + shg[t][i][l]*u(index + i)*2.0/h;
                }
            }
            eru = eru + ((solexata(xx) - uh)*(solexata(xx)-uh)*w[l]*h/2.0);
        }
        index +=nint-1; 
        erl2 = erl2 + eru;
    }
    erl2 = sqrt(erl2);
    return erl2;
}
void galerkin_continum(int nel, int nint, int nen, f h, f epslon, f gamma, contourCondition k1, contourCondition k2, f (*G)(f), f (*solexata)(f))
{
    //Inicializando a matrix e o vetor fonte
    // vector<vector<f>> K(1 + ((nint-1)*nel), vector<f>(1 + ((nint-1)*nel), 0));
    MatrixXd K = MatrixXd::Zero(1+(nint-1)*nel, 1+(nint-1)*nel);
    VectorXd F = VectorXd::Zero(1+(nint-1)*nel);
    // vector<f> F(1 + ((nint-1)*nel), 0);

    //Inicialize weight vector
    vector<f> w = gauss_weights(nint);

    //Inicializando shg
    f*** shg = create_shg(nen, nint);

    vector<f> x(nel+1, 0); 
    for(int i = 0; i <nel+1;i++)
    {
        x[i] = i*h;
    }

    for(int n = 0; n < nel; n++)
    {
        int offset = n*(nen-1); 
        for(int l=0; l<nint; l++)
        {
            f xx = 0;
            for(int i=0; i<nen; i++)
            {
                xx = xx + shg[0][i][l]*(x[n] + i*h/(nen-1));
            }

            for(int j = 0; j<nen; j++)
            {
                F(j+offset) = F(j+offset) +G(xx) * shg[0][j][l]* w[l] * h/2;

                for(int i = 0; i<nen; i++)
                {
                    K(i+offset,j+offset) = K(i+offset,j+offset)+ epslon*(shg[1][i][l] * shg[1][j][l] * 2/h * w[l]) + gamma*(shg[0][i][l]*shg[0][j][l]*w[l]*h/2);
                }
            }
        }

    }

    //Condições de contorno:
    //Dirichilet:
    if(k1.type == DIRICHLET)
    {
        F(0) = k1.value; // Condição de contorno de Dirichilet aqui
        for(int i = 1; i<nint; i++)
        {   
            F(i) -= K(i,0) * k1.value; 
            K(0, i) = 0;
            K(i, 0) = 0;
        }
        K(0,0) = 1;
    }
    else if(k1.type == NEUMANN)
    {
        F(0) -= k1.value;
    }
    
    if(k2.type == DIRICHLET)
    {
        F((nint-1)*nel) = k2.value; 

        //Os valores do lado eu não preciso fazer nada, mas com os valores abaixo é necessário fazer loucuras
        for(int i = 1; i<nint; i++) // Tirei o zero, pra ele nn fazer loucuras demais 
        {
            F(((nint-1)*nel)-i) -= K(((nint-1)*nel), ((nint-1)*nel)-i)*k2.value;
            K(((nint-1)*nel)-i, ((nint-1)*nel)) = 0;
            K(((nint-1)*nel), ((nint-1)*nel)-i) = 0;
        }
        // F(((nint-1)*nel)-1) -= K(((nint-1)*nel), ((nint-1)*nel)-1)*k2.value;
        K(((nint-1)*nel),((nint-1)*nel)) = 1;

        
    }else if(k2.type == NEUMANN)
    {
        F((nint -1)*nel) -= k2.value;
    }
    VectorXd u = K.partialPivLu().solve(F);

   
    ofstream file("output.txt");
    if(!file)
    {
        cerr<<"Erro ao abrir arquivo de saída"<<endl;
        return;
    }

    //Escrevendo X
    vector<f> xs;
    file<<u.size()<<endl;
    for(int i =0; i<u.size();i++)
    {
        xs.push_back(i*h/(nen-1));
        file<<i*h/(nen-1);
        if(i!=u.size()-1)
            file<<",";
    }
    file<<endl;

    //Escrevendo solução exata no arquivo:
    for(int i = 0; i<u.size(); i++)
    {
        file<<solexata(i*h/(nen-1));
        if(i!=u.size()-1)
            file<<",";
    }
    file<<endl;

    //Escrevendo solução no arquivo:
    for(int i = 0; i<u.size(); i++)
    {
        file<<u(i);
        if(i!=u.size()-1)
            file<<",";
    }
    file<<endl;

    file.close();
    
    cout<<"Erro solução"<<endl;
    cout<<errul2(nel, nint, nen, u, solexata, h, xs, 0)<<endl;
}

void galerking_LS(int nel, int nint, int nen, f h, f epslon, f gamma, contourCondition k1, contourCondition k2, f (*G)(f), f (*usolexata)(f), f (*psolexata)(f) ,  f del1, f del2)
{
    //Inicialize weight vector
    vector<f> w = gauss_weights(nint);
    //Inicializando shg
    f*** shg = create_shg(nen, nint);
    

    vector<f> x(nel+1, 0); 
    for(int i = 0; i <nel+1;i++)
    {
        x[i] = i*h;
    }

    //!Não sei se esses tamanhos estão certos
    VectorXd Fe = VectorXd::Zero((1+(nint-1)*nel)*2.0);
    MatrixXd Ae = MatrixXd::Zero((1+(nint-1)*nel)*2, (1+(nint-1)*nel)*2);

    for(int n = 0; n < nel; n++)
    {
        int offset = n*(nen-1); 
        for(int l=0; l<nint; l++)
        {
            f xx = 0;
            for(int i=0; i<nen; i++)
            {
                xx = xx + shg[0][i][l]*(x[n] + i*h/(nen-1));
            }

            for(int j = 0; j<nen; j++)
            {
                //G
                Fe(j+offset) += del2*G(xx)*shg[1][j][l]*2.0/h * w[l] *h/2.0; 
                //F
                Fe(j+offset + 1+(nint-1)*nel) -= G(xx) * shg[0][j][l]*w[l]*(h/2.0); 

                for(int i = 0; i<nen; i++)
                {
                    //A(linha, coluna)
                    //A 
                    Ae(i+offset, j+offset) += ((1+del1)*shg[0][i][l]*shg[0][j][l]*w[l]*h/2.0) + (del2*shg[1][i][l]*2.0/h *shg[1][j][l]*2.0/h *w[l] *h/2.0);
                    
                    //B
                    Ae(i+offset, j+offset+(1+(nint-1)*nel)) += (-1)*shg[0][j][l]*shg[1][i][l]*2.0/h*w[l]*h/2.0 + del1*shg[1][j][l]*2.0/h *shg[0][i][l]*w[l]*h/2.0;
                    //Bt                    
                    Ae(j+offset+(1+(nint-1)*nel), i+offset) += (-1)*shg[0][j][l]*shg[1][i][l]*2.0/h*w[l]*h/2.0 + del1*shg[1][j][l]*2.0/h *shg[0][i][l]*w[l]*h/2.0;
                    //C
                    Ae(i+offset + (1+(nint-1)*nel), j+offset + (1+(nint-1)*nel)) += (del1*shg[1][i][l]*(2.0/h)*shg[1][j][l]*(2.0/h))*w[l]*h/2.0;
                }
        
            }
        }
    }
    
    printmatrix(Ae);
    VectorXd u = Ae.partialPivLu().solve(Fe);

   
    ofstream file("output.txt");
    if(!file)
    {
        cerr<<"Erro ao abrir arquivo de saída"<<endl;
        return;
    }

    //Escrevendo X
    vector<f> xs;
    file<<u.size()<<endl;
    for(int i =0; i<u.size()/2;i++)
    {
        xs.push_back(i*h/(nen-1));
        file<<i*h/(nen-1);
        if(i!=(u.size()/2)-1)
            file<<",";
    }
    file<<endl;

    //Escrevendo solução exata no arquivo:
    for(int i = 0; i<u.size()/2; i++)
    {
        file<<usolexata(i*h/(nen-1));
        if(i!=(u.size()/2)-1)
            file<<",";
    }
    file<<endl;

    //Escrevendo solução de p no arquivo:
    for(int i = 0; i<(u.size()/2); i++)
    {
        file<<u(i);
        if(i!=(u.size()/2)-1)
            file<<",";
    }
    file<<endl;

    //Escrevendo a solução de u o arquivo
    for(int i = int(u.size()/2); i<u.size(); i++)
    {
        file<<u(i);
        if(i!=u.size()-1)
            file<<",";
    }
    file<<endl;

    file.close();
    
    cout<<"Erro solução u:"<<endl;
    cout<<errul2(nel, nint, nen, u.segment(0, int(u.size()/2)), usolexata, h, xs, 0)<<endl;
    
    cout<<"Erro solução p:"<<endl;
    cout<<errul2(nel, nint, nen, u.segment(int(u.size()/2),int(u.size()/2)), psolexata, h, xs, 0)<<endl;
}

void Nitsche(int nel, int nint, int nen, f h, f epslon, f gamma, f alpha, f beta, contourCondition k1, contourCondition k2, f (*G)(f), f (*solexata)(f))
{
    //Inicializando a matrix e o vetor fonte
    // vector<vector<f>> K(1 + ((nint-1)*nel), vector<f>(1 + ((nint-1)*nel), 0));
    MatrixXd K = MatrixXd::Zero(1+(nint-1)*nel, 1+(nint-1)*nel);
    VectorXd F = VectorXd::Zero(1+(nint-1)*nel);
    // vector<f> F(1 + ((nint-1)*nel), 0);

    //Inicialize weight vector
    vector<f> w = gauss_weights(nint);

    //Inicializando shg
    f*** shg = create_shg(nen, nint);

    vector<f> x(nel+1, 0); 
    for(int i = 0; i <nel+1;i++)
    {
        x[i] = i*h;
    }

    for(int n = 0; n < nel; n++)
    {
        int offset = n*(nen-1); 
        for(int l=0; l<nint; l++)
        {
            f xx = 0;
            for(int i=0; i<nen; i++)
            {
                xx = xx + shg[0][i][l]*(x[n] + i*h/(nen-1));
            }

            for(int j = 0; j<nen; j++)
            {
                F(j+offset) = F(j+offset) +G(xx) * shg[0][j][l]* w[l] * h/2;

                for(int i = 0; i<nen; i++)
                {
                    K(i+offset,j+offset) = K(i+offset,j+offset)+ epslon*(shg[1][i][l] * shg[1][j][l] * 2/h * w[l]) + gamma*(shg[0][i][l]*shg[0][j][l]*w[l]*h/2);
                }
            }
        }
    }
    F(0) += alpha*(- shg[1][0][0]*k1.value*(2.0/h)) + beta*(- shg[0][0][0]*k1.value); //Aqui eu coloco já o valor de g(a)? Ou o valor tinha que ser de v(a)? 
    K(0, 0) += shg[1][0][0]*(2.0/h)*shg[0][0][0] - (alpha*shg[1][0][0]*(2.0/h)*k1.value) - (beta*shg[0][0][0]*shg[0][0][0]);
    F((nint-1)*nel) += alpha*(shg[1][(nen-1)][nint-1]*k2.value*(2.0/h) ) + beta*(shg[0][(nen-1)][nint-1]*k2.value);
    K((nint-1)*nel, (nint-1)*nel) += (-shg[1][(nen-1)][nint-1]*(2.0/h)*shg[0][(nen-1)][nint-1]) + (alpha*shg[1][(nen-1)][nint-1]*(2.0/h)*k2.value) + (beta*shg[0][(nen-1)][nint-1]*shg[0][(nen-1)][nint-1]);

    VectorXd u = K.partialPivLu().solve(F);

   
    ofstream file("output.txt");
    if(!file)
    {
        cerr<<"Erro ao abrir arquivo de saída"<<endl;
        return;
    }

    //Escrevendo X
    vector<f> xs;
    file<<u.size()<<endl;
    for(int i =0; i<u.size();i++)
    {
        xs.push_back(i*h/(nen-1));
        file<<i*h/(nen-1);
        if(i!=u.size()-1)
            file<<",";
    }
    file<<endl;

    //Escrevendo solução exata no arquivo:
    for(int i = 0; i<u.size(); i++)
    {
        file<<solexata(i*h/(nen-1));
        if(i!=u.size()-1)
            file<<",";
    }
    file<<endl;

    //Escrevendo solução no arquivo:
    for(int i = 0; i<u.size(); i++)
    {
        file<<u(i);
        if(i!=u.size()-1)
            file<<",";
    }
    file<<endl;

    file.close();
    
    cout<<"Erro solução"<<endl;
    cout<<errul2(nel, nint, nen, u, solexata, h, xs, 0)<<endl;
}

void discontinuous_lagrange_exact(int nel, int nint, int nen, f h, f epslon, f alpha, f beta, contourCondition k1, contourCondition k2, f (*G)(f), f (*solexata)(f))
{
     //Inicializando a matrix e o vetor fonte
    // vector<vector<f>> K(1 + ((nint-1)*nel), vector<f>(1 + ((nint-1)*nel), 0));
    
    // vector<f> F(1 + ((nint-1)*nel), 0);

    //Inicialize weight vector
    vector<f> w = gauss_weights(nint);

    //Inicializando shg
    f*** shg = create_shg(nen, nint);

    VectorXd u = VectorXd::Zero(nel*2);

    vector<f> x(nel+1, 0); 
    for(int i = 0; i <nel+1;i++)
    {
        x[i] = i*h;
    }

    int index = 0;
    for(int n=0; n<nel; n++)
    {
        MatrixXd K = MatrixXd::Zero(nint, nint);
        VectorXd F = VectorXd::Zero(nint);

        int offset = n*(nen-1); 
        for(int l=0; l<nint; l++)
        {
            f xx = 0;
            for(int i =0; i<nen; i++)
            {
                xx += shg[0][i][l]*(x[n] + i*h/(nen-1));
                //Troquei o j da expressão do F por i, para eliminar um loop
                F(i) += G(xx) *shg[0][i][l]*w[l]*(h/2.0);
            }

            for(int j =0; j<nen;j++)
            {
                for(int i =0;i<nen;i++)
                {
                    K(i,j) = K(i,j)+ (shg[1][i][l] * shg[1][j][l] * (2.0/h) * w[l]);
                }
            }
        }
        for(int j=0; j<nen; j++)
        {
            F(j) += alpha*(shg[1][j][1]*(2.0/h)*solexata(x[n+1]) - shg[1][j][0]*(2.0/h)*solexata(x[n])) + beta*(shg[0][j][1]*solexata(x[n+1]) - shg[0][j][0]*solexata(x[n]));
            for(int i =0;i<nen; i++)
            {
                K(i,j) += -(shg[1][i][1]*(2.0/h)*shg[0][j][1]-shg[1][i][0]*(2.0/h)*shg[0][j][0]) + 
                alpha*(shg[1][j][1]*(2.0/h)*shg[0][i][1] - shg[1][j][0]*(2.0/h)*shg[0][i][0]) + beta*(shg[0][i][1]*shg[0][j][1] - shg[0][i][0]*shg[0][j][0]);
            }
        }

        VectorXd ue = VectorXd::Zero(nint);
        ue = K.inverse() * F;
        for(int i = 0; i<nen; i++)
        {
            u(index) =  ue(i);
            index++;
        }
    }

    for(int i=0; i<u.size(); i++)
    {
        cout<<u(i)<<",";
    }
    cout<<endl;
}

void hibrid_fem(int nel, int nint, int nen, f h, f epslon, f alpha, f beta, contourCondition k1, contourCondition k2, f (*G)(f), f (*solexata)(f))
{
     //Inicializando a matrix e o vetor fonte
    // vector<vector<f>> K(1 + ((nint-1)*nel), vector<f>(1 + ((nint-1)*nel), 0));
    
    // vector<f> F(1 + ((nint-1)*nel), 0);

    //Inicialize weight vector
    vector<f> w = gauss_weights(nint);

    //Inicializando shg
    f*** shg = create_shg(nen, nint);

    VectorXd u = VectorXd::Zero(nel*2);

    vector<f> x(nel+1, 0); 
    for(int i = 0; i <nel+1;i++)
    {
        x[i] = i*h;
    }

    int index = 0;
    for(int n=0; n<nel; n++)
    {
        MatrixXd A = MatrixXd::Zero(nint, nint);
        MatrixXd B = MatrixXd::Zero(nint, 2);
        MatrixXd D = MatrixXd::Zero(2, nint);
        VectorXd F = VectorXd::Zero(nint + 2); //Tem esse +1 por causa da equação de lambda

        for(int l=0; l<nint-1; l++)
        {
            f xx = 0;
            for(int i =0; i<nen; i++)
            {
                xx += shg[0][i][l]*(x[n] + i*h/(nen-1));
                //Troquei o j da expressão do F por i, para eliminar um loop
                F(i) += G(xx) *shg[0][i][l]*w[l]*(h/2.0);

                D(0, i) = -shg[1][i][l] + beta*shg[0][i][l];
                D(1, i) = shg[1][i][l+1] - beta*shg[0][i][l+1];

            }

            for(int j =0; j<nen;j++)
            {
                //aqui eu tenho que verificar dnv, se o primeiro item é a linha ou a coluna
                B(j, 0) = alpha*shg[1][j][l] + beta*shg[0][j][l];
                B(j, 1) = - alpha*shg[1][j][l+1] - beta*shg[0][j][l+1];

                for(int i =0;i<nen;i++)
                {
                    A(i,j) = A(i,j)+ (shg[1][i][l] * shg[1][j][l] * (2.0/h) * w[l]);
                    A(i,j) = -(shg[1][i][l+1]*(2.0/h)*shg[0][j][l+1]-shg[1][i][l]*(2.0/h)*shg[0][j][l]) + 
                    alpha*(shg[1][j][l+1]*(2.0/h)*shg[0][i][l+1] - shg[1][j][l]*(2.0/h)*shg[0][i][l]) + beta*(shg[0][i][l+1]*shg[0][j][l+1] - shg[0][i][l]*shg[0][j][l]);
                }
            }
        }
        for(int j=0; j<nen; j++)
        {
            F(j) += alpha*(shg[1][j][1]*(2.0/h)*solexata(x[n+1]) - shg[1][j][0]*(2.0/h)*solexata(x[n])) + beta*(shg[0][j][1]*solexata(x[n+1]) - shg[0][j][0]*solexata(x[n]));
        }
        Matrix2d C;
        C << -beta, 0.0, 0.0, beta;
        
        
        /* VectorXd ue = VectorXd::Zero(nint);
        ue = K.inverse() * F;
        for(int i = 0; i<nen; i++)
        {
            u(index) =  ue(i);
            index++;
        } */
    }
}

void questao1()
{   
    cout<<"Questão 1 da primeira lista de ANMEF"<<endl;
    
    for(int i = 1; i<=3; i++)
    {
        f a = 0, b = 1.5;
        int nel = pow(4, i);

        cout<<endl<<nel<<" elementos"<<endl;
        f h = (b-a)/nel;

        int nint = 5;
        int nen = nint;

        galerkin_continum(nel, nint, nen, h, 1, 0, create_contourCondition(M_PI, NEUMANN), create_contourCondition(sin(M_PI*1.5), DIRICHLET), pi2sinPiX, sinpix);
    }
    

}

void questao2()
{
    cout<<"Questão 2 da primeira lista de ANMEF"<<endl;
    f a = 0, b = 1;
    int nel = 40;
    f h = (b-a)/nel;

    int nint = 2;
    int nen = nint;
    f epslon = EPSLON;
    galerkin_continum(nel, nint, nen, h, epslon, 1, create_contourCondition(0, DIRICHLET), create_contourCondition(0, DIRICHLET), constante_1, solexata2q);
}

void lista2()
{  
    for(int i = 5; i <= 5; i++)
    {
        cout<<"Questão 1 da segunda lista de ANMEF"<<endl;
        f a = 0, b = 1;
        int nel = pow(2, i);
        f h = (b-a)/nel;
        cout<<nel<<endl;

        int nint = 2;
        int nen = nint;
        galerking_LS(nel, nint, nen, h, 0, 0, create_contourCondition(0, DIRICHLET), create_contourCondition(0, DIRICHLET), pi2sinPiX, mpicospix , sinpix, -1.0/2.0, -1.0/2.0);
    }
}

void lista3questaoA()
{
    cout<<"Questão A da terceira lista de ANMEF"<<endl;
    f a = 0, b = 1;
    int nel = 8;
    f h = (b-a)/nel;

    int nint = 2;
    int nen = nint;
    f alpha = -1;
    f beta = 1000;
    Nitsche(nel, nint, nen, h, 1, 0, alpha, beta, create_contourCondition(1, DIRICHLET), create_contourCondition(-1, DIRICHLET), pi2cospix, cospix);
}

void lista3questaoC()
{
    cout<<"Questão C da terceira lista de ANMEF"<<endl;
    f a = 0, b = 1;
    int nel = 8;
    f h = (b-a)/nel;

    int nint = 2;
    int nen = nint;
    f alpha = 1;
    f beta = 1;
    discontinuous_lagrange_exact(nel, nint, nen, h, 1, alpha, beta, create_contourCondition(1, DIRICHLET), create_contourCondition(-1, DIRICHLET), pi2cospix, cospix);
}

void lista3questaoE(){
    cout<<"Questão E da terceira lista de ANMEF"<<endl;
    f a = 0, b = 1;
    int nel = 8;
    f h = (b-a)/nel;

    int nint = 3;
    int nen = nint;
    f alpha = 1;
    f beta = 1;
    hibrid_fem(nel, nint, nen, h, 1, alpha, beta, create_contourCondition(1, DIRICHLET), create_contourCondition(-1, DIRICHLET), pi2cospix, cospix);
}


int main(){
    //Lista 1:
    //teste1();
    //questao1();
    //questao2();     
    
    //Lista 2:
    //lista2(); 

    //Lista 3:
    //lista3questaoA(); 
    //lista3questaoC(); 
    lista3questaoE();

    return 0;
}
