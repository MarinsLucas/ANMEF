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
    ofstream file("matrix.txt");
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
        pt.push_back(-1/sqrt(3.0));
        pt.push_back(1/sqrt(3.0));
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

    for(int n =0; n<nel; n++)
    {
        f eru = 0;
        for(int l = 0; l < nint; l++)
        {
            f uh = 0;
            f xx = 0;

            for(int i=0; i<nen; i++)
            {
                xx = xx + shg[t][i][l] *(x[n*(nint -1) + i]);
                    if(t==0)
                {
                    uh = uh + shg[t][i][l] * u(n*(nint -1) + i); //Tenho que selecionar o u da mesma forma do de baixo
                }else{
                    uh = uh + shg[t][i][l]*u(n*(nint -1) + i)*2.0/h;
                }
            }
            eru = eru + ((solexata(xx) - uh)*(solexata(xx)-uh)*w[l]*h/2.0);
        }
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

void galerking_LS(int nel, int nint, int nen, f h, f epslon, f gamma, contourCondition k1, contourCondition k2, f (*G)(f), f (*solexata)(f), f del2, f del1)
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
                Fe(j+offset) = 0.0; 
                //F
                Fe(j+offset + 1+(nint-1)*nel) -= G(xx) * shg[0][j][l]*w[l]*(h/2.0); 

                for(int i = 0; i<nen; i++)
                {
                    //A(linha, coluna)
                    //A 
                    Ae(i+offset, j+offset) += (shg[0][j][l] * shg[0][i][l]  + (del1 * shg[0][j][l] * shg[0][i][l])) * w[l] *h/2.0; 
                    //B
                    Ae(i+offset, j+offset+(1+(nint-1)*nel)) += (-shg[0][j][l]*shg[1][i][l] *(2.0/h) + del1*shg[1][j][l]*(2.0/h)*shg[0][i][l])*w[l]*(h/2.0);
                    //Bt                    
                    Ae(j+offset+(1+(nint-1)*nel), i+offset) += (-shg[0][j][l]*shg[1][i][l] *(2.0/h) + del1*shg[1][j][l]*(2.0/h)*shg[0][i][l])*w[l]*(h/2.0);
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
        file<<solexata(i*h/(nen-1));
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
    
    cout<<"Erro solução"<<endl;
    cout<<errul2(nel, nint, nen, u, solexata, h, xs, 0)<<endl;
}

void teste1(){
    cout<<"Teste duas condições Dirichlet do slide da Aula 1"<<endl;
    f a = 0, b = 1;
    int nel = 16; 
    f h = (b-a)/nel;
    
    int nint = 2;
    int nen = nint;
    galerkin_continum(nel, nint, nen, h, 1, 0, create_contourCondition(0, DIRICHLET), create_contourCondition(0.5, DIRICHLET), constante_1, sin);
}

void questao1()
{   
    cout<<"Questão 1 da primeira lista de ANMEF"<<endl;
    
    for(int i = 1; i<=5; i++)
    {
        f a = 0, b = 1.5;
        int nel = pow(4, i);

        cout<<endl<<nel<<" elementos"<<endl;
        f h = (b-a)/nel;

        int nint = 4;
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

void lista_questao1()
{
    cout<<"Questão 2 da primeira lista de ANMEF"<<endl;
    f a = 0, b = 1;
    int nel = 64;
    f h = (b-a)/nel;

    int nint = 2;
    int nen = nint;
    galerking_LS(nel, nint, nen, h, 0, 0, create_contourCondition(0, DIRICHLET), create_contourCondition(0, DIRICHLET), pi2sinPiX, mpicospix ,-1.0/2.0, -1.0/2.0);
}

int main(){
    //Lista 1:
    //teste1();
    //questao1();
    //questao2();     
    
    //Lista 2:
    lista_questao1(); 

    return 0;
}
