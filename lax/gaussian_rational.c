#include <stdio.h>
#include <stdlib.h>


struct gRational {
    long long rn;
    long long rd;
    long long in;
    long long id; 
};
typedef struct gRational gRational;


void stringify(gRational X, char* str) {
    sprintf(str, "%d/%d + i %d/%d\0", X.rn, X.rd, X.in, X.id);
}

void print(gRational X) {
    char str[1024];
    stringify(X, str);
    printf("%s\n", str);
}

unsigned long long lemire_gcd(unsigned long long u, unsigned long long v) {
    int shift;
    if (u == 0)
        return v;
    if (v == 0)
        return u;
    shift = __builtin_ctz(u | v);
    u >>= __builtin_ctz(u);
    do {
        unsigned long long m;
        v >>= __builtin_ctz(v);
        v -= u;
        m = (long long)v >> 63;
        u += v & m;
        v = (v + m) ^ m;
    } while (v != 0);
    return u << shift;
}

gRational red(gRational X) {
    if (X.rd == 0 || X.id == 0) {
        return X;
    }
    if (X.rn == 0) {
        X.rd = 1;

    }
    else {
        long long rgcd = lemire_gcd(abs(X.rn), abs(X.rd));
        while (rgcd != 1) {
            X.rn /= rgcd;
            X.rd /= rgcd;
            rgcd = lemire_gcd(abs(X.rn), abs(X.rd));
        }
    }
    if (X.in == 0) {
        X.id = 1;
    }
    else {
        long long igcd = lemire_gcd(abs(X.in), abs(X.id));
        while (igcd != 1) {
            X.in /= igcd;
            X.id /= igcd;
            igcd = lemire_gcd(abs(X.in), abs(X.id));
        }
    }
    return X;
}

gRational add(gRational X, gRational Y) {
    return red((gRational) {X.rn * Y.rd + X.rd * Y.rn, 
                            X.rd * Y.rd, 
                            X.in * Y.id + X.id * Y.in, 
                            X.id * Y.id});
}

gRational mul(gRational X, gRational Y) {
    return red((gRational) {X.rn * Y.rn * X.id * Y.id - X.in * Y.in * X.rd * Y.rd, 
                            X.rd * Y.rd * X.id * Y.id, 
                            X.rn * Y.in * X.id * Y.rd + X.in * Y.rn * X.rd * Y.id, 
                            X.rd * Y.id * X.id * Y.rd});
}

gRational real(gRational X) {
    return (gRational) {X.rn, 
                        X.rd, 
                        0, 
                        0};
}

gRational imag(gRational X) {
    return (gRational) {X.in, 
                        X.id, 
                        0, 
                        0};
}

gRational sqabs(gRational X) {
    return red((gRational) {X.rn * X.rn * X.id * X.id + X.rd * X.rd * X.in * X.in, 
                        X.rd * X.rd * X.id * X.id, 
                        0, 
                        0});
}

gRational inv(gRational X) {
    return red((gRational) {X.rn * (X.rd * X.rd * X.id * X.id), 
                        X.rd * (X.rn * X.rn * X.id * X.id + X.rd * X.rd * X.in * X.in),  
                        - X.in * (X.rd * X.rd * X.id * X.id), 
                        X.id * (X.rn * X.rn * X.id * X.id + X.rd * X.rd * X.in * X.in)});
}

int main() {
    gRational X = (gRational) {0, 1, 1, 2};
    gRational Y = (gRational) {3, 4, 75, 103};

    char str[1024];

    print(red((gRational) {0, 1, 1, 1}));
    print(X);
    print(Y);
    print(inv(Y));
    return 0;
}