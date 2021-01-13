#include <iostream>
#include <random>

using namespace std;

class BallsInBins
{
protected:
    ///protocol parameters and time
    int n,c,lambdan,t;
    ///randomness
    random_device rd;
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
    BallsInBins(int n, int c, int lambdan);
    ~BallsInBins();
    bool update();
    int size();
};

BallsInBins::BallsInBins(int n, int c, int lambdan) : n(n), c(c), lambdan(lambdan), gen(rd())
{
    m=vector<long long>(1,0);
    t=0;
    loads = vector<int>(c+1,0);
    loads[0] = n;
    loadHistory = vector<long long>(1000,0);
    sumLoadHistory = 0;
    takeSnapshot();
}

BallsInBins::~BallsInBins() {}

bool BallsInBins::update() {
    t++;
    deleteBalls();
    ballsGetOlder();
    createBalls();
    throwBalls();
    //bool b = updateStatistics();
    takeSnapshot();
    //return b;
    return true;
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
    vector<double> temp = vector<double>(9,0);
    //switch back 5 and 0 if needed
    temp[0] = c; // c
    temp[1] = 1; // d
    temp[2] = n; // n
    temp[3] = 1.0*lambdan/n; // lambda
    temp[4] = 0; // modified

    temp[5] = lambdan;
    for (int i=0; i<m.size(); ++i) temp[5] += (double) m[i];
    temp[5] /= n; // m/n

    for (int i=1; i <= c; ++i) temp[6] += (double) loads[i]*i;
    temp[6] /= n; // load/n

    //maximum age
    for (int i=ages.size()-1; i > 0 ; --i) {
        if (ages[i] > 0) {
            temp[7] = i;
            break;
        }
    }
    //average age
    int sum = 0;
    for (int i=0; i < ages.size(); ++i) {
        temp[8] += ages[i]*i;
        sum += ages[i];
    }
    temp[8] /= sum;

    //do whatever you want with temp here
}

int main()
{
    //start simulation
    BallsInBins simulation = BallsInBins(/*n=*/1000,/*c=*/4,/*lambda=*/0.9);
    while (simulation.update()) {}
    return 0;
}
