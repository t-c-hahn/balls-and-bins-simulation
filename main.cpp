#include <iostream>
#include <fstream>
#include <random>

using namespace std;

random_device rd();

class Snapshot
{
private:
    vector<vector<vector<double>>> data;
    int counter;
    int currentTime;
    int maxTime;
    int nStates;
public:
    ///Constructor
    Snapshot();
    ///resets the current time to 0
    void reset();
    ///access data at time step t
    vector<double> operator[](int t);
    ///add data to snapshot at current time
    ///averaging is done, time is incremented
    void append(vector<double> data);
    ///number of data points in time
    int size();
    ///current time
    int time();
};

Snapshot::Snapshot() :
    data(vector<vector<vector<double>>>()),
    counter(-1),
    currentTime(0),
    maxTime(-1),
    nStates(0) {
    //reset();
}

void Snapshot::reset() {
    data.push_back(vector<vector<double>>());
    //extend timeline to maxTime
    for (int t=currentTime; t<=maxTime; ++t)
        data[counter].push_back(vector<double>(data[counter][currentTime-1]));
    currentTime=0;
}

vector<double> Snapshot::operator[](int t) {
    vector<double> result = vector<double>(nStates,0);
    //iterate over the runs
    for (int run=0; run <= counter; ++run) {
        //if a run is too short, use the last values of the run instead
        vector<double> dataset = (run < counter || t < currentTime) ? data[run][t] : data[run][currentTime-1];
        //iterate over the states
        for (int state=0; state<dataset.size(); ++state)
            result[state]+=dataset[state]/(counter+1);
    }
    return result;
}

void Snapshot::append(vector<double> v) {
    if (nStates < v.size()) nStates = v.size();
    if (currentTime == 0) {
        counter++;
        data.push_back(vector<vector<double>>());
    }
    data[counter].push_back(vector<double>());
    data[counter][currentTime].insert(data[counter][currentTime].end(),v.begin(),v.end());
    currentTime++;
    if (currentTime>maxTime) {
        maxTime = currentTime;
        for (int run=0; run<counter; ++run)
            data[run].push_back(vector<double>(data[run][currentTime-1]));
    }
}

int Snapshot::size() {
    return maxTime;
}

int Snapshot::time() {
    return currentTime;
}

class BallsInBins
{
protected:
    ///protocol parameters and time
    int n,c,lambdan,t;
    ///randomness
    mt19937_64 gen;
    uniform_int_distribution<int> distr;

    ///balls sorted by age
    vector<long long> m;
    ///number of bins sorted by load
    vector<int> loads;

    ///age of balls that got allocated + their slot
    vector<int> ages;

    ///to determine whether the simulation has converged
    vector<long long> loadHistory;
    long long sumLoadHistory;

    ///remove one ball per bin
    virtual void deleteBalls();
    ///increase the age of the balls
    virtual void ballsGetOlder();
    ///create lambda*n new balls
    virtual void createBalls();
    ///throw the balls onto the bins
    virtual void throwBalls();
    ///update statistics at the end of the round
    ///@return: does the simulation need to continue
    virtual bool updateStatistics();
    ///store all relevant data
    virtual void takeSnapshot();

public:
    Snapshot snapshot;
    BallsInBins(int n, int c, int lambdan);
    ~BallsInBins();
    bool update();
    int size();
};

BallsInBins::BallsInBins(int n, int c, int lambdan) : n(n), c(c), lambdan(lambdan), snapshot(Snapshot())
{
    m=vector<long long>(1,0);
    t=0;
    distr=uniform_int_distribution<int>(0,n-1);
    loads = vector<int>(c+1,0);
    loads[0] = n;
    loadHistory = vector<long long>(1000,0);
    sumLoadHistory = 0;
    takeSnapshot();
}

BallsInBins::~BallsInBins() {}

bool BallsInBins::update() {
    //cout << t << endl;
    t++;
    deleteBalls();
    ballsGetOlder();
    createBalls();
    //cout << "created" << endl;
    throwBalls();
    //cout << "thrown" << endl;
    //bool b = updateStatistics();
    takeSnapshot();
    //cout << "snapshot" << endl;
    //return b;
    return t<2000;
}

void BallsInBins::deleteBalls() {
    loads[0] += loads[1];
    for (int i = 1; i < c; ++i) {
        loads[i] = loads[i+1];
    }
    loads[c] = 0;
}

void BallsInBins::ballsGetOlder() {
    if (m[m.size()-1]>0) m.resize(m.size()+1);
    for (int age = m.size()-1; age > 0; --age) {
        m[age] = m[age-1];
    }
    m[0] = 0;
}

void BallsInBins::createBalls() {
    m[0] += lambdan;
}

void BallsInBins::throwBalls() {
    int r,index;
    long long i;
    ages = vector<int>(1,0);
    for (int age = m.size()-1; age >= 0; age--) {

        long long temp_m = m[age];
        for (i=0; i<temp_m; ++i)
        {
            index = 0;
            r = distr(gen);
            while (r > loads[index] || loads[index] == 0) {
                r -= loads[index];
                index++;
            }
            if (index != c) {
                loads[index]--;
                loads[index+1]++;
                m[age]--;
                if (ages.size() <= age+index) ages.resize(age+index+1);
                ages[age+index]++;
            }
        }
    }
}

bool BallsInBins::updateStatistics() {
    //update the sum of balls in the system over the last 100 rounds
    //add the balls of the current round
    int length = 1000; //length of interval for averaging
    long long temp = 0;
    for (int i = 1; i < c; ++i) temp += loads[i];
    for (int i = 1; i < m.size(); ++i) temp += m[i];
    sumLoadHistory += temp;

    //subtract the balls 100 rounds ago
    sumLoadHistory -= loadHistory[(t-1) % length];

    //is the load 100 rounds age less than the average over the last 100 rounds?
    bool b = (loadHistory[(t-1) % length]*length < sumLoadHistory || t < length);

    //update the value 100 rounds ago with the current one
    loadHistory[(t-1) % length] = temp;

    return b;
    //return (t<1000);
}

int BallsInBins::size() {
    return n;
}

//Time, Age distribution in the cloud, -1, load distribution, -1, age at deletion of allocated balls
void BallsInBins::takeSnapshot() {
    vector<double> temp = vector<double>(5,0);
    //switch back 5 and 0 if needed
    temp[0] = c; // c

    for (int i=0; i<m.size(); ++i) temp[1] += (double) m[i];
    temp[1] /= n; // m/n

    for (int i=1; i <= c; ++i) temp[2] += (double) loads[i]*i;
    temp[2] /= n; // load/n

    //maximum age
    for (int i=ages.size()-1; i > 0 ; --i) {
        if (ages[i] > 0) {
            temp[3] = i;
            break;
        }
    }
    //average age
    int sum = 0;
    for (int i=0; i < ages.size(); ++i) {
        temp[4] += ages[i]*i;
        sum += ages[i];
    }
    temp[4] /= sum;

    snapshot.append(temp);
}

vector<double> average(Snapshot snapshot) {
    int length = 1000;
    //averaging over last length rounds
    vector<double> avg = vector<double>(snapshot[snapshot.size()-length].size(),0);
    for (int t=snapshot.size()-length; t < snapshot.size(); ++t) {
        for (int i=0; i<snapshot[t].size(); ++i) {
            avg[i] += snapshot[t][i]/length;
        }
    }
    return avg;
}

//likely bound for m/n
double boundm1(int c, double lambda) {
    return 1.0 + log(1.0/(1.0-lambda))/c;
}

//theoretical bound for m/n
double boundm2(int c, double lambda) {
    return 6.0*c + 4.0*log(1.0/(1.0-lambda))/c;
}

//theoretical bound
double boundt(int n, int c, double lambda) {
    return log2(log2(n)) + c + log(1.0/(1.0-lambda))/c;
}

void run1() {
    ofstream file("plot1.csv");
    if (!file.is_open()) return;

    vector<BallsInBins*> sims = vector<BallsInBins*>(10);
    sims[0] = new BallsInBins(/*n=2^15*/32768,/*c=*/1,/*lambda*n=3/4*n*/24576);
    sims[1] = new BallsInBins(32768,2,24576);
    sims[2] = new BallsInBins(32768,3,24576);
    sims[3] = new BallsInBins(32768,4,24576);
    sims[4] = new BallsInBins(32768,5,24576);
    sims[5] = new BallsInBins(32768,1,32736);
    sims[6] = new BallsInBins(32768,2,32736);
    sims[7] = new BallsInBins(32768,3,32736);
    sims[8] = new BallsInBins(32768,4,32736);
    sims[9] = new BallsInBins(32768,5,32736);
    
    vector<vector<double>> avgs = vector<vector<double>>(10);
    for (int i=0; i<10;++i) {
        while (sims[i]->update()) {}
        avgs[i] = average(sims[i]->snapshot);
    }

    //Concatenate all data
    file << "m(c);lambda=3/4;ln(4)/c+1;4ln(4)+6c;lambda=1023/1024;ln(1024)/c+1;4ln(1024)+6c" << endl; //header
    for (int i = 0; i < 5; ++i) {
        file << to_string(i+1) << ";"
             << to_string(avgs[i][1]) << ";"
             << to_string(boundm1(i+1,0.75)) << ";"
             << to_string(boundm2(i+1,0.75)) << ";"
             << to_string(avgs[5+i][1]) << ";"
             << to_string(boundm1(i+1,1023.0/1024.0)) << ";"
             << to_string(boundm2(i+1,1023.0/1024.0)) << endl;
    }
}

void run2() {
    ofstream file("plot2.csv");
    if (!file.is_open()) return;

    vector<BallsInBins*> sims = vector<BallsInBins*>(20);
    int slack = 32768;
    for (int i = 0; i < 10; ++i) {
        slack /= 2;
        sims[i] = new BallsInBins(/*n=2^15*/32768,/*c=*/1,/*lambda*n=1/2*n*/32768-slack);
        sims[i+10] = new BallsInBins(/*n=2^15*/32768,/*c=*/3,/*lambda*n=1/2*n*/32768-slack);
    }
    
    vector<vector<double>> avgs = vector<vector<double>>(20);
    for (int i=0; i<20;++i) {
        while (sims[i]->update()) {}
        avgs[i] = average(sims[i]->snapshot);
    }

    //Concatenate all data
    file << "c=1;bound(c=1);c=3;bound(c=3)" << endl; //header
    double x = 1;
    for (int i = 0; i < 10; ++i) {
        x *= 2;
        file << to_string(avgs[i][1]) << ";"
             << to_string(boundm1(1,1-1/x)) << ";"
             << to_string(avgs[10+i][1]) << ";"
             << to_string(boundm1(3,1.0-1.0/x)) << endl;
    }
}

void run3() {
    ofstream file("plot3.csv");
    if (!file.is_open()) return;

    vector<BallsInBins*> sims = vector<BallsInBins*>(15);
    sims[0] = new BallsInBins(/*n=2^15*/32768,/*c=*/1,/*lambda*n=3/4*n*/24576);
    sims[1] = new BallsInBins(32768,2,24576);
    sims[2] = new BallsInBins(32768,3,24576);
    sims[3] = new BallsInBins(32768,4,24576);
    sims[4] = new BallsInBins(32768,5,24576);
    sims[5] = new BallsInBins(32768,1,32736);
    sims[6] = new BallsInBins(32768,2,32736);
    sims[7] = new BallsInBins(32768,3,32736);
    sims[8] = new BallsInBins(32768,4,32736);
    sims[9] = new BallsInBins(32768,5,32736);
    sims[10] = new BallsInBins(32768,1,32764);
    sims[11] = new BallsInBins(32768,2,32764);
    sims[12] = new BallsInBins(32768,3,32764);
    sims[13] = new BallsInBins(32768,4,32764);
    sims[14] = new BallsInBins(32768,5,32764);    

    
    vector<vector<double>> avgs = vector<vector<double>>(15);
    for (int i=0; i<15;++i) {
        while (sims[i]->update()) {}
        avgs[i] = average(sims[i]->snapshot);
    }

    //Concatenate all data
    file << "tau(c);max(3/4);avg(3/4);max(1023/1024);avg(1023/1024);avg(8191/8192);max(8191/8192);bound(4);bound(1024);bound(8192)" << endl; //header
    for (int i = 0; i < 5; ++i) {
        file << to_string(i+1) << ";"
             << to_string(avgs[i][3]) << ";"
             << to_string(avgs[i][4]) << ";"
             << to_string(avgs[5+i][3]) << ";"
             << to_string(avgs[5+i][4]) << ";"
             << to_string(avgs[10+i][3]) << ";"
             << to_string(avgs[10+i][4]) << ";"
             << to_string(boundt(32768,i+1,0.75)) << ";"
             << to_string(boundt(32768,i+1,1023.0/1024.0)) << ";"
             << to_string(boundt(32768,i+1,8191.0/8192.0)) << endl;
    }
}

void run4() {
    ofstream file("plot4.csv");
    if (!file.is_open()) return;

    vector<BallsInBins*> sims = vector<BallsInBins*>(20);
    int slack = 32768;
    for (int i = 0; i < 10; ++i) {
        slack /= 2;
        sims[i] = new BallsInBins(/*n=2^15*/32768,/*c=*/1,/*lambda*n=1/2*n*/32768-slack);
        sims[i+10] = new BallsInBins(/*n=2^15*/32768,/*c=*/3,/*lambda*n=1/2*n*/32768-slack);
    }
    
    vector<vector<double>> avgs = vector<vector<double>>(20);
    for (int i=0; i<20;++i) {
        while (sims[i]->update()) {}
        avgs[i] = average(sims[i]->snapshot);
    }

    //Concatenate all data
    file << "lambda;max(c=1);avg(c=1);max(c=3);avg(c=3);bound(c=1);bound(c=3)" << endl; //header
    double x = 1;
    for (int i = 0; i < 10; ++i) {
        x *= 2;
        file << to_string(1.0-1.0/x) << ";"
             << to_string(avgs[i][3]) << ";"
             << to_string(avgs[i][4]) << ";"
             << to_string(avgs[10+i][3]) << ";"
             << to_string(avgs[10+i][4]) << ";"
             << to_string(boundt(32768,1,1.0-1.0/x)) << ";"
             << to_string(boundt(32768,3,1.0-1.0/x)) << endl;
    }
}

int main()
{
    //start simulations for figure 1,2,3 and 4
    run1();
    run2();
    run3();
    run4();
    return 0;
}

