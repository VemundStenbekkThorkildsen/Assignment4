#include <iostream>
#include <armadillo>
#include <random>

using namespace std;
using namespace arma;

random_device rd;
mt19937 randomEngine(rd());
uniform_real_distribution<double> uniformDist(0.0,1.0);

ofstream outFile, outFile2, outFile3, outFile4, outFile5, outFile6;

double random(){
    return uniformDist(randomEngine);
}

void toFile(double Mtemp, double Etemp, double T, int acceptance){
    outFile << Mtemp << ", " << endl;
    outFile2 << Etemp << ", " << endl;
    outFile4 << acceptance << ", " << endl;
}

void toFile2(int acceptance, double T){
    outFile3 << T << ", " << acceptance << endl;
}

void toFile3(int cycles, double E, double M){
    outFile5 << cycles << "," << E << endl;
    outFile6 << cycles << "," << M << endl;
}
mat randomMatrix(mat &A, int L){

    for(int i = 0; i < L; i++){
        for(int j = 0; j < L; j++){
            if(random() <= 0.5){
                A(i,j) = 1;
            }
            else{
                A(i,j) = -1;
            }
        }
    }
    return A;
}


void MP(int L, mat &A, vec prob, int &acceptance, double &E, double &Mtemp, double &E_2, double &Etemp, double &M, double &M_2){
    int xp, xn, yp, yn;
    double deltaE;
    for(int x = 0; x < L; x++){
        for(int y = 0; y < L; y++){
            int xran=round(random()*(L-1));
            int yran=round(random()*(L-1));

            xp = (xran + 1) % L;
            yp = (yran + 1) % L;
            xn = (xran - 1 + L) % L;
            yn = (yran - 1 + L) % L;

            deltaE = 2.0 * A(xran,yran) * (  A(xp,yran) + A(xn,yran) + A(xran, yn) + A(xran, yp));

            if(random() <= prob(deltaE + 8)){
                A(xran, yran) *= -1;
                Etemp += deltaE;

                Mtemp += 2*A(xran,yran);
                acceptance += 1;

            }
        }
    }

}


void openFiles(double T){
        string outFileName = "magnetization" + to_string(T) + ".txt";
        string outFileName2 = "energy" + to_string(T) + ".txt";
        string outFileName4 = "acceptanceMCup" + to_string(T) + ".txt";
        outFile.open(outFileName);
        outFile2.open(outFileName2);
        outFile4.open(outFileName4);
}

void openFiles2(){
    string outFileName3 = "TempAcceptanceUp.txt";
    outFile3.open(outFileName3);
}

void openFiles3(double T){
    string outFileName5 = "CumEnergyRAND" + to_string(T) +".txt";
    string outFileName6 = "CumMagneRAND" + to_string(T) +".txt";
    outFile5.open(outFileName5);
    outFile6.open(outFileName6);
}
int main()
{
    //openFiles2();
    //double T=1.0;
    for(double T = 2.4; T <= 2.4; T+=0.6){
    double k = 1.0; //J/K
    double beta = 1.0/(k*T);
    double J = 1.0;

    int mcs = 1000000;

    int L = 20;

    openFiles(T);
    openFiles3(T);

    mat A = ones(L,L);
    //randomMatrix(A, L);



    //Initial values of the temporary energy
    double Etemp = 0;
    for(int x = 0; x<L; x++){
        for(int y = 0; y<L;y++){

            int xn = (x - 1 + L) % L;
            int yn = (y - 1 + L) % L;

            Etemp -= A(x,y) * (  A(xn,y) + A(x, yn));
        }
    }

    double E = 0;
    double E_2 = 0;
    double Mtemp = accu(A);
    double M = 0;
    double M_2= 0;

    vec prob(17);
    for(int i=-8; i <= 8; i+=4){
        prob(i+8) = 0;
    }
    for(int i=-8; i <= 8; i+=4){
        prob(i + 8) = exp(-i/T);
    }

    int acceptance = 0;
    for(int cycles=0;cycles<=mcs;cycles++){
        MP(L, A, prob, acceptance, E, Mtemp, E_2, Etemp, M, M_2);


       if(cycles >= 100000){
            toFile(Mtemp, Etemp, T, acceptance);
            toFile3(cycles,E,M);
        }
            E += Etemp;
            E_2 += Etemp*Etemp;
            M += abs(Mtemp);
            M_2 += Mtemp*Mtemp;




    }

    //toFile2(acceptance, T);

    // Mean energy
    double averegeE = E/(mcs);

    // Mean magnetism
    double averegeM = (M)/(mcs);

    // Specific heat
    double averegeESquared = E_2/(mcs+1);
    double heat = (beta*(averegeESquared - (averegeE*averegeE)))/T;

    // Susceptibility
    double averegeMSquared = (M_2/(mcs+1));
    double sus = beta*(averegeMSquared - averegeM*averegeM);



    cout << endl << "Average energy: " << averegeE << " and the square of average energy: " << averegeESquared << " while T = " << to_string(T) << endl;
    cout << "Average magnetization: " << averegeM <<  " while T = " << to_string(T) << endl;
    cout << "Specific heat: " << heat << " while T = " << to_string(T) << endl;
    cout << "Susceptibility: " << sus << " while T = " << to_string(T) << endl;
    cout << "Variance Energy:" << (heat*T*T)/(L*L) << "while T=" << to_string(T) << endl;

    outFile.close();
    outFile2.close();
    outFile4.close();
    outFile5.close();
    outFile6.close();
    }
    //outFile3.close();



    return 0;
}
