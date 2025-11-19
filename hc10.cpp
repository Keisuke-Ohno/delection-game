#include <stdio.h>
#include <stdlib.h>
#include <bitset>
#include <iostream>
#include <map>
#include <cstdint>

#define N 10/*vertexの数*/

/*可変ビット幅での1bit左循環シフト*/
uint32_t rotl1(uint32_t x, unsigned int width) {
    width = width % 32; // 念のため32ビット以上の値を正規化
    return ((x << 1) | (x >> (width - 1))) & ((1u << width) - 1);
}

int size = 1 << N; /*配列checkの大きさ*/


/*bitsetの大小を定義する関数(10bitまで)*/


struct Comparer {
    bool operator()(const std::bitset<1024> &b1, const std::bitset<1024> &b2) const {
        uint64_t part1[16] = {0}, part2[16] = {0};

        // 64ビットずつ整数に詰める
        for (int block = 0; block < 16; ++block) {
            for (int i = 0; i < 64; ++i) {
                int bit_index = block * 64 + i;
                if (b1[bit_index]) part1[block] |= (1ULL << i);
                if (b2[bit_index]) part2[block] |= (1ULL << i);
            }
        }

        // 上位ブロックから比較（辞書順）
        for (int block = 15; block >= 0; --block) {
            if (part1[block] != part2[block])
                return part1[block] < part2[block];
        }

        return false; // 全て同じ → a == b
    }
};





int calc_total = 0; /*Calcを実行した数を格納する*/

std::map<std::bitset<1024>, int, Comparer> m; 

void Remove(int *Check, int num); /*Checkから取り除く*/
void Print_Check(int *Check); /*Checkを出力*/
int Continue(int *Check, int size); /*全て取り除かれたか判定*/
int* Create_Check(void); /*Checkを作成*/
void Print_bits(int a); /*要素を出力*/
int count_bits(int n); /*引数を2進数にしたときに幾つ１が出てくるか*/
int Calc(int *Check); /*再帰関数でg-valueを計算*/
int mex(int *arr, int size); /*引数の配列の要素の中で最大の非負整数を返す*/
int p_uniform(int *Check);

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

    uint32_t pp = 15;
    
    for(i = 0; i < N; i++) {
        Check[pp] = 1;
        pp = rotl1(pp, N);
    }
    //Remove(Check,1);
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
    std::bitset<1024> b_check;

/*現在のCheckの内容をbitsetにする*/
    for(i = 0; i < size; i++) {
        if(Check[i] == 1) {
            b_check.set(i);
        } else {
            b_check.reset(i);
        }
    }

/*parity uniformを満たしていない局面を表示*/
/*頂点しかないような局面の場合は出力させないための処理*/
    for(i = 0; i < size; i++) {
        if(count_bits(i) > 1 && Check[i] == 1) {
            flag = 1;
            break;
        }
    }
    if(p_uniform(Check)== 0 && flag == 1) {
        printf("\n");
        printf("@@@@@@@@@\n");
        Print_Check(Check);
        printf("@@@@@@@@@\n");
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


    /*式を使って求めた結果と計算した結果が一致しない局面を表示*/
    /*
    if(parity_uni != g_value) {
        printf("\n");
        printf("************************\n");
        printf("V = %d, E = %d\n", v, e);
        printf("g_value = %d, parity_uniform = %d\n", g_value, parity_uni);
        Print_Check(Check);
        printf("************************");
        printf("\n");
    }*/

    

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

/*その局面がparity uniformがどうかを求める*/
int p_uniform(int *Check)
{
    int i,j;
    int OK = 1;

    for(i = 0; i < size; i++) {
        for(j = 0; j < size; j++) {
            if(Check[j] == 1 && count_bits(j) > 1) {
                if((count_bits(i & j) % 2) == 0) {
                    OK = 0;
                    break;
                }
            }        
        }
        if(OK == 1) {
            return i;
        }
        OK = 1;
    }
    return 0;
}

    
