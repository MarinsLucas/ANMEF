/* TO COMPILE:
 g++ *.cpp *.h -o ./output.exe
 
 */
#include<vector>
#include<cmath>
#include<iostream>
#include<eigen3/Eigen/Dense>

#define f float

#define DIRICHLET 0
#define NEUMANN 1


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

f G(f x)
{
    return 1;
}

void printvector(vector<f> v)
{
    for(int i = 0; i<v.size(); i++)
    {
        cout<<v[i]<<",";
    }
    cout<<endl;
}

void printmatrix(vector<vector<f>> m)
{
    for(int i = 0; i<m.size(); i++)
    {  
        cout<<"[";
        for (int j = 0; j<m[0].size(); j++)
        {
            cout<<m[i][j]<<",";
        }
        cout<<"]"<<endl;
    }
}

// Função para realizar a decomposição LU
void decomposeLU(std::vector<std::vector<f>>& A, std::vector<std::vector<f>>& L, std::vector<std::vector<f>>& U) {
    int n = A.size();

    // Inicializa L e U com zeros
    L = std::vector<std::vector<f>>(n, std::vector<f>(n, 0));
    U = std::vector<std::vector<f>>(n, std::vector<f>(n, 0));

    // Preenche U com os elementos de A
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            U[i][j] = A[i][j];
        }
    }

    // Preenche a diagonal de L com 1s
    for (int i = 0; i < n; i++) {
        L[i][i] = 1;
    }

    // Realiza a decomposição LU
    for (int k = 0; k < n - 1; k++) {
        for (int i = k + 1; i < n; i++) {
            L[i][k] = U[i][k] / U[k][k];
            for (int j = k; j < n; j++) {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }
    }
}

// Função para resolver um sistema linear Ax = b usando a decomposição LU
std::vector<f> solveUsingLU(std::vector<std::vector<f>>& A, std::vector<f>& b) {
    int n = A.size();
    std::vector<std::vector<f>> L, U;

    // Realiza a decomposição LU
    decomposeLU(A, L, U);

    // Resolve Ly = b
    std::vector<f> y(n);
    for (int i = 0; i < n; i++) {
        f sum = 0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = (b[i] - sum) / L[i][i];
    }

    // Resolve Ux = y
    std::vector<f> x(n);
    for (int i = n - 1; i >= 0; i--) {
        f sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }

    return x;
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

    vector<f> pt(nint, 0);
    switch (nint)
    {
    case 2:
        pt.push_back(-1/sqrt(3));
        pt.push_back(1/sqrt(3));
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
        
        default:
            break;
        }
    }

    return shg; 
}

void fem(int nel, int nint, int nen, f h, f epslon, f gamma, contourCondition k1, contourCondition k2)
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
        int offset = n + (nen-2); 
        for(int l=0; l<nint; l++)
        {
            f xx = 0;
            for(int i=0; i<nen; i++)
            {
                xx = xx + shg[0][i][l]*(x[n] + i*h/nen);
            }

            for(int j = 0; j<nen; j++)
            {
                F(j+offset) = F(j+offset) + G(xx) * shg[0][j][l]* w[l] * h/2;

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
        F(1) -= K(1,0) * k1.value; 
        K(0, 1) = 0;
        K(1, 0) = 0;
        K(0,0) = 1;
    }
    else if(k1.type == NEUMANN)
    {
        F(0) += k1.value;
    }
    
    if(k2.type == DIRICHLET)
    {
        F((nint-1)*nel) = k2.value; 
        F(((nint-1)*nel)-1) -= K(((nint-1)*nel), ((nint-1)*nel)-1)*k2.value;
        K(((nint-1)*nel),((nint-1)*nel)) = 1;
        K(((nint-1)*nel)-1, ((nint-1)*nel)) = 0;
        K(((nint-1)*nel), ((nint-1)*nel)-1) = 0;
    }else if(k2.type == NEUMANN)
    {
        F((nint -1)*nel) += k2.value;
    }

    VectorXd u = K.partialPivLu().solve(F);
    cout<<"Solução"<<endl<<u<<endl;
}

void teste1(){
    cout<<"Teste duas condições Dirichlet do slide da Aula 1"<<endl;
    f a = 0, b = 1;
    int nel = 4; 
    f h = (b-a)/nel;
    
    int nint = 2;
    int nen = nint;
    
    fem(nel, nint, nen, h, 1, 0, create_contourCondition(0, DIRICHLET), create_contourCondition(0.5, DIRICHLET));
}

void questao1()
{
    cout<<"Questão 1 da primeira lista de ANMEF"<<endl;
    f a = 0, b = 1.5;
    int nel = 16;
    f h = (b-a)/nel;

    int nint = 2;
    int nen = nint;

    fem(nel, nint, nen, h, 1, 0, create_contourCondition(M_PI, DIRICHLET), create_contourCondition(sin(M_PI*1.5), NEUMANN));

}

void questao2()
{
    cout<<"Questão 2 da primeira lista de ANMEF"<<endl;
    f a = 0, b = 1;
    int nel = 4;
    f h = (b-a)/nel;

    int nint = 2;
    int nen = nint;
    f epslon = 0.001;
    fem(nel, nint, nen, h, epslon, 1, create_contourCondition(0, DIRICHLET), create_contourCondition(0, DIRICHLET));
}

int main(){
    teste1();
    questao2();

    return 0;
} 
