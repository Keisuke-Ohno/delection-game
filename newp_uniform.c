#include <stdio.h>
#include <stdlib.h>
#define N 4 /*vertexの数*/


int size = 1 << N; /*配列checkの大きさ*/
int gauss_gf2(int E, int V, int A[E][N+1], int x[N]);
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
    Remove(Check,1);
    Remove(Check,2);
    Remove(Check,4);
    Remove(Check,8);

    
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
   int tmp,answer;
   

    for(i = 0; i < size; i++){
        if(Check[i] == 1 && count_bits(i) == 1) {
            V++;
        } else if(Check[i] == 1 && count_bits(i) > 1) {
            Edge[E] = i;
            E++;

        }
    }
    
   int A[E][N+1];
   int x[V];

   printf("V = %d, E = %d\n", V, E);

   for(i = 0; i < E; i++) {
    printf("%d ", Edge[i]);
   }
   printf("\n");



    for(i = 0; i < E; i++) {
        tmp = 1;
        for(j = 0; j < N; j++) {
            k = Edge[i] & tmp;
            if(k > 0) { A[i][j] = 1; }
            else { A[i][j] = 0; }
            tmp = tmp << 1; 
        }
    }

    for(i = 0; i < E; i++) {
        A[i][N] = 1;
    }

    

    printf("[");
    for(i = 0; i < E; i++) {
        
        printf("[");
        for(j = 0; j < N+1; j++) {
            printf("%d, ", A[i][j]);
        }
        printf("]");
        printf("\n");
    }
    printf("]\n");

    answer = gauss_gf2(E, N, A, x);

    if(answer == 0) {
        printf("answer = {");
        for(i = 0; i < N; i++) {
            printf("%d, ", x[i]);
        }
        printf("}\n");
    } else {
        printf("解なし\n");
    }
    


}

int gauss_gf2(int E, int V, int A[E][V+1], int x[V])
{
    int row = 0;

    for (int col = 0; col < V && row < E; col++) {

        // ピボット探索
        int pivot = -1;
        for (int i = row; i < E; i++) {
            if (A[i][col]) {
                pivot = i;
                break;
            }
        }
        if (pivot == -1) continue;

        // 行交換
        if (pivot != row) {
            for (int j = col; j <= V; j++) {
                int tmp = A[row][j];
                A[row][j] = A[pivot][j];
                A[pivot][j] = tmp;
            }
        }

        // 他の行を掃き出す（GF(2)なので XOR）
        for (int i = 0; i < E; i++) {
            if (i != row && A[i][col]) {
                for (int j = col; j <= V; j++) {
                    A[i][j] ^= A[row][j];
                }
            }
        }

        row++;
    }

    // 解のチェックと取り出し
    for (int i = 0; i < V; i++) x[i] = 0;

    for (int i = 0; i < E; i++) {
        int lead = -1;
        for (int j = 0; j < V; j++) {
            if (A[i][j]) {
                lead = j;
                break;
            }
        }
        if (lead == -1) {
            if (A[i][V]) return -1; // 矛盾 → 解なし
        } else {
            x[lead] = A[i][V];
        }
    }

    return 0; // 解あり（自由変数は0にしている）
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



