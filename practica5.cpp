#include <iostream>
#include <vector>
#include <string>
#include <omp.h>
#include <chrono>
#include <cstring>

#define N 100000000 // Tamaño de la secuencia
const std::string PATTERN = "ATGC"; // Patrón a buscar
const int P_LEN = PATTERN.length();
using namespace std;

// Función para simular una secuencia grande
void generate_sequence(vector<char>& seq) {
    // Rellenamos con bases aleatorias para simular una carga de trabajo real
    // Insertamos el patrón cerca del final para forzar la búsqueda larga
    srand(time(0));
    for (int i = 0; i < N - P_LEN; ++i) {
        seq[i] = "ATGC"[rand() % 4];
    }
    // Aseguramos que el patrón exista al final para que la búsqueda no falle
    for (int i = 0; i < P_LEN; ++i) {
        seq[N - 10000 + i] = PATTERN[i]; 
    }
}

void run_search(int num_threads, const char* schedule_type, int chunk_size) {
    vector<char> dna_sequence(N);
    generate_sequence(dna_sequence);
    long long first_index = -1; 
    omp_set_num_threads(num_threads);

    auto start = chrono::high_resolution_clock::now();

    // ================================
    //     Selección del Schedule
    // ================================
    if (strcmp(schedule_type, "Secuencial") == 0) {

        for (int i = 0; i < N - P_LEN; ++i) {

            if (first_index != -1) continue;

            bool match = true;
            for (int j = 0; j < P_LEN; ++j)
                if (dna_sequence[i + j] != PATTERN[j]) { match = false; break; }

            if (match) {
                if (first_index == -1 || i < first_index)
                    first_index = i;
            }
        }

    } else if (strcmp(schedule_type, "Static") == 0) {

        #pragma omp parallel for schedule(static, chunk_size)
        for (int i = 0; i < N - P_LEN; ++i) {
            if (first_index != -1) continue;
            bool match = true;
            for (int j = 0; j < P_LEN; ++j)
                if (dna_sequence[i + j] != PATTERN[j]) { match = false; break; }

            if (match) {
                #pragma omp critical
                if (first_index == -1 || i < first_index)
                    first_index = i;
            }
        }

    } else if (strcmp(schedule_type, "Dynamic") == 0) {

        #pragma omp parallel for schedule(dynamic, chunk_size)
        for (int i = 0; i < N - P_LEN; ++i) {
            if (first_index != -1) continue;
            bool match = true;
            for (int j = 0; j < P_LEN; ++j)
                if (dna_sequence[i + j] != PATTERN[j]) { match = false; break; }

            if (match) {
                #pragma omp critical
                if (first_index == -1 || i < first_index)
                    first_index = i;
            }
        }

    } else if (strcmp(schedule_type, "Guided") == 0) {

        #pragma omp parallel for schedule(guided, chunk_size)
        for (int i = 0; i < N - P_LEN; ++i) {
            if (first_index != -1) continue;
            bool match = true;
            for (int j = 0; j < P_LEN; ++j)
                if (dna_sequence[i + j] != PATTERN[j]) { match = false; break; }

            if (match) {
                #pragma omp critical
                if (first_index == -1 || i < first_index)
                    first_index = i;
            }
        }

    } else if (strcmp(schedule_type, "Auto") == 0) {

        #pragma omp parallel for schedule(auto)
        for (int i = 0; i < N - P_LEN; ++i) {
            if (first_index != -1) continue;
            bool match = true;
            for (int j = 0; j < P_LEN; ++j)
                if (dna_sequence[i + j] != PATTERN[j]) { match = false; break; }

            if (match) {
                #pragma omp critical
                if (first_index == -1 || i < first_index)
                    first_index = i;
            }
        }

    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;

    cout << "Hilos: " << num_threads 
         << ", Schedule: " << schedule_type 
         << " (Chunk: " << chunk_size << ")"
         << ", Tiempo: " << elapsed.count() << " s" 
         << ", Posición: " << first_index << endl;
}

int main() {
    cout << "--- Búsqueda de Patron en ADN (" << PATTERN << ") ---" << endl;
    cout << "---------------------------------------------------" << endl;

    int num_hilos = omp_get_max_threads(); //Numero de Hilos de la VM

    // 1. Ejecución Secuencial (1 Hilo) 
    cout << "--- 1. Ejecución Secuencial ---" << endl; run_search(1, "Secuencial", 0); 
    // 2. Ejecución Paralela - Dynamic con Chunk Pequeño 
    cout << "\n--- 2. Ejecución Dynamic ---" << endl; run_search(num_hilos, "Dynamic", 10); 
    // 3. Ejecución Paralela - Dynamic con Chunk Grande 
    cout << "\n--- 3. Ejecución Dynamic ---" << endl; run_search(num_hilos, "Dynamic", N); 
    // 4. Ejecución Paralela - Static con Chunk Pequeño 
    cout << "\n--- 4. Ejecución Static ---" << endl; run_search(num_hilos, "Static", 10); 
    // 5. Ejecución Paralela - Static con Chunk Grande 
    cout << "\n--- 5. Ejecución Static ---" << endl; run_search(num_hilos, "Static", N); 
    // 6. Ejecución Paralela - Guided con Chunk Pequeño 
    cout << "\n--- 6. Ejecución Guided ---" << endl; run_search(num_hilos, "Guided", 10); 
    // 7. Ejecución Paralela - Guided con Chunk Grande 
    cout << "\n--- 7. Ejecución Guided ---" << endl; run_search(num_hilos, "Guided", N);
    // 8. Ejecución Paralela - Auto
    cout << "\n--- 8. Ejecución Auto ---" << endl; run_search(num_hilos, "Auto", 0); 
    
    return 0;
}
