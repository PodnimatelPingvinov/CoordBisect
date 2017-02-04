#include <stdio.h>

int main(int argc, char **argv)
{
    if(argc != 2){
        printf("Usage: ./view <file>\n");
        return 1;
    }
      
    FILE *input = fopen(argv[1], "rb");
    if(input == NULL){
        printf("File '%s' can't be opened.\n", argv[1]);
        return 1;
    }
    
    int n1, n2, i, j, domain;
    float x, y;
    fread(&n1, sizeof(n1), 1, input);
    fread(&n2, sizeof(n2), 1, input);
    
    for(int k = 0; k < n1 * n2; k++){
        fread(&i, sizeof(i), 1, input);
        fread(&j, sizeof(j), 1, input);
        fread(&x, sizeof(x), 1, input);
        fread(&y, sizeof(y), 1, input);
        fread(&domain, sizeof(domain), 1, input);
        printf("%d %d %f %f %d\n", i, j, x, y, domain);
    }
    
    fclose(input);
    
    return 0;
}
