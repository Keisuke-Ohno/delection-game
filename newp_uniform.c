#include <stdio.h>
#include <stdlib.h>
#define N 6 /*vertexの数*/

int size = 1 << N; /*配列checkの大きさ*/

uint32_t rotl1(uint32_t x, unsigned int width) {
    width = width % 32; // 念のため32ビット以上の値を正規化
    return ((x << 1) | (x >> (width - 1))) & ((1u << width) - 1);
}

void Remove(int *Check, int num);

void Print_bits(int a)
{
    int i;
    int b;
    printf(" {");
    for(i = 0; i < N; i++) {
        b = (a>>i) & 1;
        if(b == 1) {
            printf("%d, ", i+1);
        }

    }
    printf("}");
}
void Print_Check(int *Check)
{
    int i;

    printf("\n");
    for(i = 0; i < size; i++) {
        if(Check[i] == 1) {
            printf("%d = ", i);
            Print_bits(i);
            printf("\n");
        }
    }
}
int count_bits(int n)
{
  unsigned int x = (unsigned int)n;
    int sum = 0;

    while (x) {
        x &= (x - 1);  // 最下位の1を消す
        sum++;
    }
    return sum;
}
int* Create_Check(void)
{
    int i;
    
     int* Check = (int*)malloc(size * sizeof(int));  // 動的に配列を確保
    if (Check == NULL) {
        printf("メモリの確保に失敗しました\n");
        return NULL;
    }
    /*配列の初期化 */
    for(i = 0; i < size; i++) {
        Check[i] = 0;
    }

    for(i = 0; i < size; i++) {
        if(count_bits(i) == 1) {
           Check[i] = 1;
        }
    }
    uint32_t pp = 15;
    
    for(i = 0; i < N; i++) {
        Check[pp] = 1;
        pp = rotl1(pp, N);
    }
    //Remove(Check,387);
    //Remove(Check,30);
    //Remove(Check,57);
    return Check;
} 
void Remove(int *Check, int num) 
{
    int i;
    int tmp;

    Check[num] = 0;
    
    for(i = 1; i < size; i++) {
        tmp = i & num;
        if(tmp == num) {
            Check[i] = 0;
        }
    }
}

void p_uniform(int *Check)
{
   int i, j, k;
   int V=0,E=0;
   int Edge[N];
   int tmp;
   

    for(i = 0; i < size; i++){
        if(Check[i] == 1 && count_bits(i) == 1) {
            V++;
        } else if(Check[i] == 1 && count_bits(i) > 1) {
            Edge[E] = i;
            E++;

        }
    }
    
   int A[E][V];

   printf("V = %d, E = %d\n", V, E);

   for(i = 0; i < E; i++) {
    printf("%d ", Edge[i]);
   }
   printf("\n");



    for(i = 0; i < E; i++) {
        tmp = 1;
        for(j = 0; j < V; j++) {
            k = Edge[i] & tmp;
            if(k > 0) { A[i][j] = 1; }
            else { A[i][j] = 0; }
            tmp = tmp << 1; 
        }
    }

    for(i = 0; i < E; i++) {
        for(j = 0; j < V; j++) {
            printf("%d, ", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");


}

int main(void)
{
    int* Check = Create_Check();

    Print_Check(Check);

    p_uniform(Check);


    printf("size = %d\n", size);

    

    
    

    free(Check);
    return 0;
}



