int bfpf(int N,int a){
    int exp = a;
    for(int i = 1; i < N; ++i){
        if(exp % N == 1){
            return i;
            break;
        }
        exp *= a;
		exp = exp % N;
    }
    return -1;
}
int main(){
	printf("period is:%d", bfpf(1007,529));
}
