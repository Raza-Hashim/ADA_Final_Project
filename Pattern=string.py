import time
import numpy as np
import matplotlib.pyplot as plt
import random
import string

NO_OF_CHARS = 256

def badCharHeuristic(string, size):
    badChar = [-1] * NO_OF_CHARS
    for i in range(size):
        badChar[ord(string[i])] = i
    return badChar

def preprocess_strong_suffix(shift, bpos, pat, m):
    i = m
    j = m + 1
    bpos[i] = j
    while i > 0:
        while j <= m and pat[i - 1] != pat[j - 1]:
            if shift[j] == 0:
                shift[j] = j - i
            j = bpos[j]
        i -= 1
        j -= 1
        bpos[i] = j

def preprocess_case2(shift, bpos, pat, m):
    j = bpos[0]
    for i in range(m + 1):
        if shift[i] == 0:
            shift[i] = j
        if i == j:
            j = bpos[j]

def BM(text, pat, n, m):
    badChar = badCharHeuristic(pat, m)
    bpos = [0] * (m + 1)
    shift = [0] * (m + 1)
    preprocess_strong_suffix(shift, bpos, pat, m)
    preprocess_case2(shift, bpos, pat, m)
    s = 0
    while s <= n - m:
        j = m - 1
        while j >= 0 and pat[j] == text[s + j]:
            j -= 1
        if j < 0:
            print("Pattern occurs at shift = {}".format(s))
            s += shift[0]
        else:
            bad_char_shift = j - badChar[ord(text[s + j])]
            good_suffix_shift = shift[j + 1]
            s += max(1, max(bad_char_shift, good_suffix_shift))

def naive_string_matcher(T, P, n, m):
    matches = []
    for s in range(n - m + 1):
        if T[s:s + m] == P:
            matches.append(s)
    return matches


def compute_prefix_function(pattern, m):
    prefix = [0] * m
    k = 0
    for q in range(1, m):
        while k > 0 and pattern[k] != pattern[q]:
            k = prefix[k-1]
        if pattern[k] == pattern[q]:
            k += 1
        prefix[q] = k
    return prefix

def kmp_matcher(text, pattern, n, m):
    prefix = compute_prefix_function(pattern, m)
    q = 0
    for i in range(n):
        while q > 0 and pattern[q] != text[i]:
            q = prefix[q-1]
        if pattern[q] == text[i]:
            q += 1
        if q == m:
            pass  # Do nothing, remove print statement for timing
            q = prefix[q-1]

def RABIN_KARP_MATCHER(text, pattern, n, m, q=389, d=256):
    h = pow(d, m - 1, q)
    p = 0
    t = 0
    for i in range(m):
        p = (d * p + ord(pattern[i])) % q
        t = (d * t + ord(text[i])) % q

    matches = []
    for s in range(n - m + 1):
        if p == t and pattern == text[s:s + m]:
            matches.append(s)

        if s < n - m:
            t = (d * (t - ord(text[s]) * h) + ord(text[s + m])) % q
            if t < 0:
                t += q

    return matches




def time_naive_string_matcher(texts, patterns, num_trials=10):
    times = []
    for text, pattern in zip(texts, patterns):
        avg_time = 0
        for _ in range(num_trials):
            start_time = time.time()
            naive_string_matcher(text, pattern, len(text), len(pattern))
            end_time = time.time()
            avg_time += (end_time - start_time)
        avg_time /= num_trials
        times.append(avg_time)
    return times

def time_kmp_matcher(texts, patterns, num_trials=10):
    times = []
    for text, pattern in zip(texts, patterns):
        avg_time = 0
        for _ in range(num_trials):
            start_time = time.time()
            kmp_matcher(text, pattern, len(text), len(pattern))
            end_time = time.time()
            avg_time += (end_time - start_time)
        avg_time /= num_trials
        times.append(avg_time)
    return times

def time_rk_string_matcher(texts, patterns, num_trials=10):
    times = []
    for text, pattern in zip(texts, patterns):
        avg_time = 0
        for _ in range(num_trials):
            start_time = time.time()
            RABIN_KARP_MATCHER(text, pattern, len(text), len(pattern))
            end_time = time.time()
            avg_time += (end_time - start_time)
        avg_time /= num_trials
        times.append(avg_time)
    return times

def time_bm_matcher(texts, patterns, num_trials=10):
    times = []
    for text, pattern in zip(texts, patterns):
        avg_time = 0
        for _ in range(num_trials):
            start_time = time.time()
            BM(text, pattern, len(text), len(pattern))
            end_time = time.time()
            avg_time += (end_time - start_time)
        avg_time /= num_trials
        times.append(avg_time)
    return times

# Generate random strings and patterns
num_strings = 10
# strings, patterns, string_lengths = generate_random_strings(num_strings)

strings=[]
patterns=[]
string_lengths=[]

string_="ABCD"

for i in range(num_strings):
    l =  random.randint(2500, 250000)
    s = string_*l
    strings.append(s)
    patterns.append(s)
    string_lengths.append(len(s))

# Time taken by naive string matcher
naive_times = time_naive_string_matcher(strings, patterns)

# Time taken by KMP matcher
kmp_times = time_kmp_matcher(strings, patterns)

# Time taken by Rabin-Karp matcher
rk_times = time_rk_string_matcher(strings, patterns)

# Time taken by Boyer-Moore matcher
bm_times = time_bm_matcher(strings, patterns)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(string_lengths, naive_times, marker='o', label='Naive Matcher')
plt.plot(string_lengths, kmp_times, marker='x', label='KMP Matcher')
plt.plot(string_lengths, rk_times, marker='x', label='RK Matcher')
plt.plot(string_lengths, bm_times, marker='x', label='BM Matcher')
plt.title('Time Taken when Pattern=String')
plt.xlabel('Text Length')
plt.ylabel('Time (seconds)')
plt.legend()
plt.grid(True)
plt.show()
