#include<stdio.h>
#include<gmp.h>
#include <windows.h>
#include <Wincrypt.h>

//b504f333f9de6484597d89b3754abe9f1d6f60ba893ba84ced17ac85833399154afc83043ab8a2c3a8b1fe6fdc83db390f74a85e439c7b4a780487363dfa2768      更号2
#define root2512 "b504f333f9de6484597d89b3754abe9f1d6f60ba893ba84ced17ac85833399154afc83043ab8a2c3a8b1fe6fdc83db390f74a85e439c7b4a780487363dfa2768"
//4AFB0CCC06219B7BA682764C8AB54160E2909F4576C457B312E8537A7CCC66EAB5037CFBC5475D3C574E0190237C24C6F08B57A1BC6384B587FB78C9C205D898      number
#define numberuproot2512 "4AFB0CCC06219B7BA682764C8AB54160E2909F4576C457B312E8537A7CCC66EAB5037CFBC5475D3C574E0190237C24C6F08B57A1BC6384B587FB78C9C205D898"

int generatepq ( mpz_t* p,mpz_t* q, const mpz_t e);
int WinRand(unsigned char* a,unsigned int lenth);
int generatep( mpz_t p, const mpz_t p1, const mpz_t p2, const mpz_t e);
int lucustest(const mpz_t N);

int main()
{
    mpz_t  p, q ,e;
    mpz_init(p);
    mpz_init(q);
    mpz_init(e);
    mpz_set_str(e,"65537",10);
    generatepq(&p,&q,e);
    getchar();
    return 0;
}

int generatepq ( mpz_t* p_ptr,mpz_t* q_ptr, const mpz_t e)
{
    mpz_t xp1, xp2, xq1, xq2 ,p1 ,p2 ,q1 ,q2;
    mpz_init(xp1);
    mpz_init(xp2);
    mpz_init(xq1);
    mpz_init(xq2);
    mpz_init(p1);
    mpz_init(p2);
    mpz_init(q1);
    mpz_init(q2);
    char rand1[26];
    WinRand(rand1,25);
    rand1[25]='\0';
    mpz_set_str(xp1, rand1, 16);            //16--hex
    gmp_printf("xp1 = %Zx\n", xp1);
    mpz_nextprime(p1,xp1);
        /*
    while(mpz_probab_prime_p(xp1,27) == 0)          //miller robin
    {
        mpz_add_ui(xp1,xp1,1);
        gmp_printf("xp1 = %Zx\n", xp1);
    }
    gmp_printf("xp1 = %Zx\n", xp1);         */    
    gmp_printf("p1 = %Zx\n", p1);
    WinRand(rand1,25);
    rand1[25]='\0';
    mpz_set_str(xp2, rand1, 16); 
    gmp_printf("xp2 = %Zx\n", xp2);
    mpz_nextprime(p2,xp2);
    gmp_printf("p2 = %Zx\n", p2);
    //generate p
    generatep(*p_ptr,p1,p2,e);

}

int WinRand(unsigned char* a,unsigned int lenth)
{
    //unsigned char a[1024];
    HCRYPTPROV Rnd;
    LPCSTR UserName = "MyKeyContainer";
    if(CryptAcquireContext(&Rnd, UserName, NULL, PROV_RSA_FULL, 0)) {
        //printf("A cryptographic context with the %s key container ",
        //UserName);
        //printf("has been acquired.\n\n");
    } else {
        if (GetLastError() == NTE_BAD_KEYSET) {
            if(CryptAcquireContext(
                &Rnd,
                UserName,
                NULL,
                PROV_RSA_FULL,
                CRYPT_NEWKEYSET)) {
                //printf("A new key container has been created.\n");
            } else {
                printf("Could not create a new key container.\n");
                //exit(1);
                return 0;
            }
        } else {
            printf("A cryptographic service handle could not be "
            "acquired.\n");
            //exit(1);
            return 0;
        }
    }
 
    if (CryptGenRandom(Rnd, lenth, (BYTE*)(a)))
    {
         if (CryptReleaseContext(Rnd,0)) {
          //printf("The handle has been released.\n\n");
            //return 1;
         } else {
            printf("The handle could not be released.\n\n");
            //return 0;
        }
        int i;
        char temp;
        a[0] = a[0]|(0x80);                 //random top is 1xxx xxxx
        for(i=0;i<lenth;i++)
        {  
            sprintf(&temp,"%02x",a[i]);
            a[i] = temp;
        }  
        return 1;
    } 
    else
    {
        if (CryptReleaseContext(Rnd,0)) {
          //printf("The handle has been released.\n\n");
            //return 1;
        } else {
            printf("The handle could not be released.\n\n");
            //return 0;
        }
        printf("Generate random number failed.\n\n");
        return 0;
    }
 
   
}

int generatep(mpz_t p, const mpz_t p1, const mpz_t p2, const mpz_t e)
{
    //random xp
    unsigned char * randxp;
    randxp = (unsigned char *)malloc(129);
    randxp[128] = '\0';
    mpz_t xp, * const_root2_ptr;
    const_root2_ptr = (mpz_t*)malloc(sizeof(mpz_t));
    //const_numup_ptr = (mpz_t*)malloc(sizeof(mpz_t));
    mpz_init(xp);
    mpz_init(*const_root2_ptr);
    //mpz_init(*const_numup_ptr);
    mpz_set_str(* const_root2_ptr,root2512,16);
    //mpz_set_str(* const_numup_ptr,numberuproot2512,16);
    do
    {
        WinRand(randxp,128);
        mpz_set_str(xp,randxp,16);
    } while (mpz_cmp(xp,* const_root2_ptr)<0);
    gmp_printf("xp = %Zx\n",xp);
    free(randxp);
    free(const_root2_ptr);
    //free(const_numup_ptr);

    // calculate Rp
    mpz_t temp1, temp2, rp, p1p2;
    mpz_init(temp1);
    mpz_init(temp2);
    mpz_init(p1p2);
    mpz_init(rp);
    mpz_mul(p1p2, p1, p2);

    mpz_invert(temp1,p2,p1);
    mpz_invert(temp2,p1,p2);
    mpz_mul(temp1,temp1,p2);
    mpz_mul(temp2,temp2,p1);
    mpz_sub(rp, temp1, temp2);

    if(mpz_cmp_ui(rp,0)<0)
        mpz_add(rp,rp,p1p2);
    gmp_printf("rp = %Zx\n", rp);
    mpz_clear(temp2);

    //calculate y0
    mpz_t yi;
    mpz_init(yi);
    mpz_mod(temp1,xp,p1p2);     //xp mod p1p2
    mpz_sub(temp1,rp,temp1);    //rp - xp mod p1p2
    mpz_add(yi,xp,temp1);       //yi = xp + (rp - xp mod p1p2)
    if(mpz_cmp(yi,xp)<0)
        mpz_add(yi, yi, p1p2);
    gmp_printf("yi = %Zx\n", yi);

    //yi + j * p1p2
    if(mpz_odd_p(e))
    {
        mpz_t i;
        mpz_init(i);
        mpz_set_ui(i,0);
        do
        {
            mpz_mul(temp1,i,p1p2);
            mpz_add(yi, yi , temp1);
        } while ((mpz_probab_prime_p(yi,8) == 0));
        
    }
    else
    {
        ;
    }

}

int lucustest(const mpz_t N)
{
    mpz_t D, P ,Q, temp1;
    mpz_init(D);
    mpz_init(P);
    mpz_init(Q);
    mpz_init(temp1);
    mpz_set_si(D,5);
    while (mpz_jacobi(D,N)!=-1)
    {
        mpz_add_ui(D,D,2);
        mpz_neg(D,D);
    }
    gmp_printf("lucas D = %Zx\n",D);
    mpz_set_si(P,1);
    mpz_ui_sub(Q, 1, D);
    mpz_div_ui(Q , Q, 4);

    //begin the steps
    mpz_mul_ui(temp1, Q, 4);
    mpz_ui_sub(D, 1, temp1);
    
    
}