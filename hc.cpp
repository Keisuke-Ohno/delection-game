#include <stdio.h>
#include <stdlib.h>
#include <bitset>
#include <iostream>
#include <map>
#include <cstdint>

#define N 8/*vertexの数*/
#define bset_size 4096/*bitsetの大きさ(2^N)*/

/*可変ビット幅での1bit左循環シフト*/
uint32_t rotl1(uint32_t x, unsigned int width) {
    width = width % 32; // 念のため32ビット以上の値を正規化
    return ((x << 1) | (x >> (width - 1))) & ((1u << width) - 1);
}

int size = 1 << N; /*配列checkの大きさ*/


/*bitsetの大小を定義する関数*/

/* 今回は12bitをサポートするため std::bitset<4096> を使います。
   4096 = 2^12, 64ブロックの 64-bit ワードへ分割して辞書順比較します。 */
struct Comparer {
    bool operator()(const std::bitset<bset_size> &b1, const std::bitset<bset_size> &b2) const {
        const int BLOCKS = bset_size / 64; // 64
        uint64_t part1[BLOCKS];
        uint64_t part2[BLOCKS];

        // 初期化
        for (int block = 0; block < BLOCKS; ++block) {
            part1[block] = 0;
            part2[block] = 0;
        }

        // 64ビットずつ整数に詰める
        for (int block = 0; block < BLOCKS; ++block) {
            for (int i = 0; i < 64; ++i) {
                int bit_index = block * 64 + i;
                if (b1[bit_index]) part1[block] |= (1ULL << i);
                if (b2[bit_index]) part2[block] |= (1ULL << i);
            }
        }

        // 上位ブロックから比較（辞書順）
        for (int block = BLOCKS - 1; block >= 0; --block) {
            if (part1[block] != part2[block])
                return part1[block] < part2[block];
        }

        return false; // 全て同じ → a == b
    }
};

int calc_total = 0; /*Calcを実行した数を格納する*/

std::map<std::bitset<bset_size>, int, Comparer> m; 

void Remove(int *Check, int num); /*Checkから取り除く*/
void Print_Check(int *Check); /*Checkを出力*/
int Continue(int *Check, int size); /*全て取り除かれたか判定*/
int* Create_Check(void); /*Checkを作成*/
void Print_bits(int a); /*要素を出力*/
int count_bits(int n); /*引数を2進数にしたときに幾つ１が出てくるか*/
int Calc(int *Check); /*再帰関数でg-valueを計算*/
int mex(int *arr, int size); /*引数の配列の要素の中で最大の非負整数を返す*/
int p_uniform(int *Check); /*与えられた局面がp_unoformを満たしているかを判定*/
int gauss_gf2(int E, int V, int **A, int *x);

int main(void)
{
    int* Check = Create_Check();

    Print_Check(Check);


    printf("size = %d\n", size);

    
   printf("result = %d\n", Calc(Check)); 
   printf("calc_total = %d\n", calc_total);
   printf("mapsize=%d\n", (int)m.size());
   
    free(Check);
    return 0;
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

    //辺の大きさを決める　大きさ４ 15、大きさ６　63
    uint32_t pp = 15; 
    
    for(i = 0; i < N; i++) {
        Check[pp] = 1;
        pp = rotl1(pp, N);
    }
    //Remove(Check,15);
    //Remove(Check,30);
    //Remove(Check,60);
    //Remove(Check,64);

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

int Continue(int *Check)
{
    int i;

    for(i = 0; i < size; i++) {
        if(Check[i] == 1) {
            return 1;
        }
    }
    return 0;
}

int Calc(int *Check)
{
    calc_total++;
    int i, j, k;
    int flag = 0;
    int g_value; /*与えられた局面のg値を格納*/
    int N_P;/**/
    int Check_copy[size];/**/
    int mex_o[size];/*与えられた局面から１手進んだ局面のg値を格納する配列*/
    int v,e,parity_uni;
    std::bitset<bset_size> b_check;

/*現在のCheckの内容をbitsetにする*/
    for(i = 0; i < size; i++) {
        if(Check[i] == 1) {
            b_check.set(i);
        } else {
            b_check.reset(i);
        }
    }





/*もし前に全く同じ局面のg-valueを計算していたのならばその時の結果を返す*/
    if(m.count(b_check) != 0) {
        return m.at(b_check);
    }

    k = 0;
    for(i = 0; i < size; i++) {
        if(Check[i] == 1) {
            for(j = 0; j < size; j++) {
                Check_copy[j] = Check[j];
            }

            Remove(Check_copy, i);
            N_P = Calc(Check_copy);

            mex_o[k] = N_P;
            k++;

        }

    //mapに追加
    }
    g_value = mex(mex_o, k);
    m.insert(std::make_pair(b_check, g_value));

    v = 0;
    e = 0;
    for(i = 0; i < size; i++){
        if(count_bits(i) == 1 && Check[i] == 1) {
            v++;
        } else if(count_bits(i) > 1 && Check[i] == 1) {
            e++;
        }
    }
    v = v % 2;
    e = e % 2;
    parity_uni = v ^ (2*e);

    /*parity uniformを満たしていない局面を表示*/
/*頂点しかないような局面の場合は出力させないための処理*/
    for(i = 0; i < size; i++) {
        if(count_bits(i) > 1 && Check[i] == 1) {
            flag = 1;
            break;
        }
    }
    /* p_uniform now returns -1 on no solution, or returns a witness i (0..size-1) */
    if(p_uniform(Check) == -1 && flag == 1) {
        printf("@@@@@@@@@");
        Print_Check(Check);
        printf("@@@@@@@@@\n");
    }


    return g_value;

}

int mex(int *arr, int size) {
    for (int i = 0; ; i++) {   // 0 から順番に探す
        int found = 0;
        for (int j = 0; j < size; j++) {
            if (arr[j] == i) {
                found = 1;
                break;
            }
        }
        if (!found) {
            return i;  // 含まれていない最小の非負整数を返す
        }
    }
}

int p_uniform(int *Check)
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
    
   int **A;
   int x[N];

    A = (int **)malloc(sizeof(int *) * E);

    for (int i = 0; i < E; i++) {
        A[i] = (int *)malloc(sizeof(int) * (N + 1));
    }

    //辺と頂点の数を表示、配列edgeの中身を表示
    //printf("V = %d, E = %d\n", V, E);
    /*for(i = 0; i < E; i++) {
     printf("%d ", Edge[i]);
    }
    printf("\n");*/



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

    //隣接行列の表示
    /*printf("[");
    for(i = 0; i < E; i++) {
        
        printf("[");
        for(j = 0; j < V+1; j++) {
            printf("%d, ", A[i][j]);
        }
        printf("]");
        printf("\n");
    }
    printf("]\n");
    */

    answer = gauss_gf2(E, N, A, x);

    for (int i = 0; i < E; i++) {
    delete[] A[i];   // 各行を解放
    }
    delete[] A;          // 行ポインタ配列を解放        

    return answer;


}

int gauss_gf2(int E, int V, int **A, int *x)
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