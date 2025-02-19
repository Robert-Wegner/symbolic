#include <stdio.h>
#include <stdlib.h>

#define MAX_LENGTH 60000

struct Mon {
    char pre;   // bit 1, 2 are the sign: index in ['+', '-', 'i', 'j']
                // bit 3 is 0 for number or 1 for constant
                // bit 4 bit is empty
                // bit 5, 6 are the value,
                // bit 7, 8 are empty
    char expr[15]; // all bits 0 is 'd', otherwise first 4 bits are variable and last 4 bits are derivatives
};
typedef struct Mon Mon;

struct Mon2 {
    char pre;   // bit 1, 2 are the sign: index in ['+', '-', 'i', 'j']
                // bit 3 is 0 for number or 1 for constant
                // bit 4 bit is 0 for empty or 1 for constant
                // bit 5, 6 are the value,
                // bit 7, 8 are the other value
    char expr[28]; // all bits 0 is 'd', otherwise first 4 bits are variable and last 4 bits are derivatives
};
typedef struct Mon2 Mon2;

struct Sum {
    Mon args[1024];
};
typedef struct Sum Sum;

struct Sum2 {
    Mon args[1024];
};
typedef struct Sum2 Sum2;

void stringify_mon(Mon mon, char* str) {
    str[0] = 'h';
    str[1] = '\0';
};

int main() {
    printf("hlooe");
    return 0;
}
