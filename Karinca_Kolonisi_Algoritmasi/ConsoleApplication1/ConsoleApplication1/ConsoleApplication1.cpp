#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <random>
#include <algorithm>
#include <chrono>

using namespace std;

// ========================================================
//  Schwefel Fonksiyonu (2 boyut: x=(x0, x1))
// ========================================================
double schwefelFunction(const vector<double>& x) {
    double result = 418.9829 * x.size();  // x.size()=2
    for (auto xi : x) {
        result -= xi * sin(sqrt(fabs(xi)));
    }
    return result;
}

// ========================================================
//  Rastgele sayı üretimi için yardımcı sınıf
// ========================================================
class RandomGenerator {
public:
    RandomGenerator() {
        rng.seed(static_cast<unsigned long>(time(nullptr)));
    }

    double uniformReal(double minVal, double maxVal) {
        uniform_real_distribution<double> dist(minVal, maxVal);
        return dist(rng);
    }

private:
    mt19937 rng;
};

// ========================================================
//  Ana Program
// ========================================================
int main() {
    // Süre ölçümü başlangıcı
    auto startTime = chrono::high_resolution_clock::now();

    // Parametreler
    const int dimension = 2;
    const int numAnts = 1000;
    const int maxIterations = 1000;
    const double desiredFitness = 1e-6;

    // Feromon parametreleri
    double pheromoneEvapRate = 0.2;
    double pheromoneInit = 1.0;
    double alpha = 1.0;  // Feromonun etkisi
    double beta = 2.0;   // Heuristik bilginin etkisi

    // Çözüm aralığı
    const double minCoord = 400.0;
    const double maxCoord = 450.0;

    // Rastgele sayı üreteci
    RandomGenerator randGen;

    // Feromonlar ve heuristik bilgiler
    vector<double> pheromones(dimension, pheromoneInit);
    vector<vector<double>> distances(numAnts, vector<double>(dimension, 0.0));

    // Global en iyi çözüm bilgisi
    vector<double> bestSolution(dimension, 0.0);
    double bestCost = numeric_limits<double>::infinity();

    // Ana döngü
    for (int iter = 0; iter < maxIterations; iter++) {
        double iterationBestCost = numeric_limits<double>::infinity();
        vector<double> iterationBestSolution(dimension);

        // Karıncaların çözümleri ve maliyetleri
        vector<vector<double>> candidates(numAnts, vector<double>(dimension, 0.0));
        vector<double> costs(numAnts);

        // Karıncaların geçiş olasılıklarını hesaplama
        for (int ant = 0; ant < numAnts; ant++) {
            vector<double> candidate(dimension);

            for (int d = 0; d < dimension; d++) {
                // Feromon ve heuristik bilgi
                double tau = pheromones[d];
                double eta = 1.0 / (maxCoord - minCoord);  // Heuristik bilgi (örnek olarak sabit)

                // Geçiş olasılığı formülü
                double numerator = pow(tau, alpha) * pow(eta, beta);
                double denominator = 0.0;

                // Alternatif yollar için toplam
                for (int k = 0; k < dimension; k++) {
                    double tau_k = pheromones[k];
                    double eta_k = 1.0 / (maxCoord - minCoord);
                    denominator += pow(tau_k, alpha) * pow(eta_k, beta);
                }

                double probability = numerator / denominator;

                // Rastgele seçim
                double randomVal = randGen.uniformReal(0.0, 1.0);
                if (randomVal <= probability) {
                    // Seçilen yol
                    candidate[d] = randGen.uniformReal(minCoord, maxCoord);
                }
            }

            candidates[ant] = candidate;
            costs[ant] = schwefelFunction(candidate);

            // Iterasyonun en iyisini güncelle
            if (costs[ant] < iterationBestCost) {
                iterationBestCost = costs[ant];
                iterationBestSolution = candidate;
            }
        }

        // Global en iyiyi güncelle
        if (iterationBestCost < bestCost) {
            bestCost = iterationBestCost;
            bestSolution = iterationBestSolution;
        }

        // Feromon güncellemesi
        for (int d = 0; d < dimension; d++) {
            pheromones[d] *= (1.0 - pheromoneEvapRate);  // Buharlaşma
            pheromones[d] += 1.0 / (1.0 + iterationBestCost);  // Iterasyon en iyisi katkısı
        }

        // Ara bilgi
        if (iter % 100 == 0) {
            cout << "[Iter " << iter << "] Best Cost = " << bestCost << endl;
        }

        // Erken durdurma
        if (bestCost < desiredFitness) {
            cout << "\nErken durdurma: hedef fitness (" << desiredFitness << ") sağlandı.\n";
            break;
        }
    }

    // Süre ölçümü sonu
    auto endTime = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsedTime = endTime - startTime;

    // Sonuç
    cout << "\n** En iyi sonuc (cost): " << bestCost << "\n";
    cout << "** En iyi (x, y): ";
    for (auto val : bestSolution) {
        cout << val << " ";
    }
    cout << endl;

    // Geçen süre
    cout << "\nAlgoritmanin calismasi icin gecen sure: " << elapsedTime.count() << " saniye.\n";

    return 0;
}