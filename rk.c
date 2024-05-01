#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <vector>

using namespace std;

#define d 256

void search(const char pat[], const char txt[], int q, std::ofstream& outputFile) {
    int M = strlen(pat);
    int N = strlen(txt);
    int i, j;
    int p = 0;
    int t = 0;
    int h = 1;

    for (i = 0; i < M - 1; i++)
        h = (h * d) % q;

    for (i = 0; i < M; i++) {
        p = (d * p + pat[i]) % q;
        t = (d * t + txt[i]) % q;
    }

    for (i = 0; i <= N - M; i++) {
        if (p == t) {
            for (j = 0; j < M; j++) {
                if (txt[i + j] != pat[j])
                    break;
            }
            if (j == M)
                outputFile << "Pattern found at index " << i << std::endl;
        }
        if (i < N - M) {
            t = (d * (t - txt[i] * h) + txt[i + M]) % q;
            if (t < 0)
                t = (t + q);
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <text> <pattern>" << std::endl;
        return 1;
    }

    const char *txt = argv[1];
    const char *pat = argv[2];

    int q = 101; // You can choose any prime number for q

    std::ofstream outputFile("output.txt");
    search(pat, txt, q, outputFile);

    return 0;
}

