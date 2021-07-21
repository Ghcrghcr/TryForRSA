//Inspired by gmpy2 : https://pypi.org/project/gmpy2/
//Used the The GNU Multiple Precision Arithmetic Library : https://gmplib.org/

#include<stdio.h>
#include<gmp.h>
#include <windows.h>
#include <Wincrypt.h>

//b504f333f9de6484597d89b3754abe9f1d6f60ba893ba84ced17ac85833399154afc83043ab8a2c3a8b1fe6fdc83db390f74a85e439c7b4a780487363dfa2768   更号2
#define root2512 "b504f333f9de6484597d89b3754abe9f1d6f60ba893ba84ced17ac85833399154afc83043ab8a2c3a8b1fe6fdc83db390f74a85e439c7b4a780487363dfa2768"
//4AFB0CCC06219B7BA682764C8AB54160E2909F4576C457B312E8537A7CCC66EAB5037CFBC5475D3C574E0190237C24C6F08B57A1BC6384B587FB78C9C205D898      number
#define numberuproot2512 "4AFB0CCC06219B7BA682764C8AB54160E2909F4576C457B312E8537A7CCC66EAB5037CFBC5475D3C574E0190237C24C6F08B57A1BC6384B587FB78C9C205D898"

int generatepq ( mpz_t* p,mpz_t* q, const mpz_t e);
int WinRand(unsigned char* a,unsigned int lenth);
int generatep( mpz_t p, const mpz_t p1, const mpz_t p2, const mpz_t e, const mpz_t xp);
//int lucustest(const mpz_t N);
int GMPY_mpz_is_selfridge_prp(const mpz_t n);
int GMPY_mpz_is_lucas_prp(const mpz_t p, const mpz_t q,const mpz_t n);

int main()
{
    mpz_t  p, q ,e;
    mpz_init(p);
    mpz_init(q);
    mpz_init(e);
    mpz_set_str(e,"65537",10);
    generatepq(&p,&q,e);
    printf("Type any thing to quit\n");
    getchar();
    return 0;
}

int generatepq ( mpz_t* p_ptr,mpz_t* q_ptr, const mpz_t e)
{
    mpz_t xp1, xp2 ,p1 ,p2;
    mpz_init(xp1);
    mpz_init(xp2);
    mpz_init(p1);
    mpz_init(p2);
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

    //random xp
    unsigned char * randxp;
    randxp = (unsigned char *)malloc(129);
    randxp[128] = '\0';
    mpz_t xp, * const_root2_ptr;
    const_root2_ptr = (mpz_t*)malloc(sizeof(mpz_t));
    mpz_init(xp);
    mpz_init(*const_root2_ptr);
    mpz_set_str(* const_root2_ptr,root2512,16);
    do
    {
        WinRand(randxp,128);
        mpz_set_str(xp,randxp,16);
    } while (mpz_cmp(xp,* const_root2_ptr)<0);
    gmp_printf("xp = %Zx\n",xp);
    free(randxp);
    free(const_root2_ptr);
    //random xp end

    //generate p
    if(generatep(*p_ptr,p1,p2,e,xp))
        gmp_printf("Find p = %Zx\n\n",*p_ptr);
    else
        printf("Occurred error in finding p!\n\n");

    //generate q
    mpz_t temptest;
    mpz_init(temptest);
    do
    {
        mpz_set_str(xp1, rand1, 16);            //16--hex
        gmp_printf("xq1 = %Zx\n", xp1);
        mpz_nextprime(p1,xp1);
            /*
        while(mpz_probab_prime_p(xp1,27) == 0)          //miller robin
        {
            mpz_add_ui(xp1,xp1,1);
            gmp_printf("xp1 = %Zx\n", xp1);
        }
        gmp_printf("xp1 = %Zx\n", xp1);         */    
        gmp_printf("q1 = %Zx\n", p1);
        WinRand(rand1,25);
        rand1[25]='\0';
        mpz_set_str(xp2, rand1, 16); 
        gmp_printf("xq2 = %Zx\n", xp2);
        mpz_nextprime(p2,xp2);
        gmp_printf(" q2 = %Zx\n", p2);

        //random xq
        unsigned char * randxq;
        randxq = (unsigned char *)malloc(129);
        randxq[128] = '\0';
        mpz_t xq, * const_root2_ptr, tempxpq;
        const_root2_ptr = (mpz_t*)malloc(sizeof(mpz_t));
        mpz_init(xq);
        mpz_init(*const_root2_ptr);
        mpz_init(tempxpq);
        mpz_set_str(* const_root2_ptr,root2512,16);
        do
        {
            WinRand(randxq,128);
            mpz_set_str(xq,randxq,16);
            mpz_sub(tempxpq,xp,xq);
            mpz_abs(tempxpq,tempxpq);
            mpz_div_2exp(tempxpq,tempxpq,412);
        } while ((mpz_cmp(xq,* const_root2_ptr)<0)||(mpz_cmp_si(tempxpq,0)<=0));
        gmp_printf("xq = %Zx\n",xq);
        mpz_clear(tempxpq);
        free(randxq);
        free(const_root2_ptr);
        //random xp end

        if(generatep(*q_ptr,p1,p2,e,xq))
            gmp_printf("Find first q = %Zx\n",*q_ptr);
        else
            printf("Occurred error in finding q!\n\n");
        mpz_sub(temptest,*p_ptr,*q_ptr);
        mpz_abs(temptest,temptest);
        mpz_div_2exp(temptest,temptest,412);
    } while (mpz_cmp_si(temptest,0)<=0);
    gmp_printf("Find fitable q = %Zx\n\n",*q_ptr);
    
    return 1;
    
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

int generatep(mpz_t p, const mpz_t p1, const mpz_t p2, const mpz_t e, const mpz_t xp)
{
    

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
    gmp_printf("y0 = %Zx\n", yi);

    //yi + j * p1p2
    if(mpz_odd_p(e))
    {
        mpz_t i ,yiSubOne, gcdans;
        mpz_init(i);
        mpz_init(yiSubOne);
        mpz_init(gcdans);
        mpz_set_ui(i,0);
        do
        {
            mpz_mul(temp1,i,p1p2);
            mpz_add(yi, yi , temp1);
            mpz_sub_ui(yiSubOne,yi,1);
            mpz_gcd(gcdans,yiSubOne,e);

            mpz_add_ui(i,i,1);
        } while ((mpz_probab_prime_p(yi,8) == 0)||!(GMPY_mpz_is_selfridge_prp(yi))|| (mpz_cmp_ui(gcdans,1) != 0));
        mpz_set(p, yi);
    }
    else
    {
        return 0;
    }
    return 1;
}

/*
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
    mpz_t U,V;
    mpz_init(U);
    mpz_init(V);
    mpz_mul_ui(temp1, Q, 4);
    mpz_ui_sub(D, 1, temp1);
    mpz_set_ui(U,1);
    mpz_set(V, P);
    size_t NUM = 64*mpz_size(N);
    size_t i;
    for(i = NUM-1; i >= 0;i--)
    {
        mpz_mul()
    }
    
}
*/

int GMPY_mpz_is_lucas_prp(const mpz_t p, const mpz_t q,const mpz_t n)
{
    mpz_t zD, res, index, uh, vl ,vh, ql, qh, tmp;
    mpz_init(zD);
    mpz_init(res);
    mpz_init(index);
    mpz_init(uh);
    mpz_init(vl);
    mpz_init(vh);
    mpz_init(ql);
    mpz_init(qh);
    mpz_init(tmp);
    int result;

    /* Check if p*p - 4*q == 0. */
    mpz_mul(zD, p, p);
    mpz_mul_ui(tmp, q, 4);
    mpz_sub(zD, zD, tmp);
    if (mpz_sgn(zD) == 0) {
        printf("invalid values for p,q in is_lucas_prp()");
        goto cleanup;
    }

    if (mpz_cmp_ui(n, 2) < 0) {
        result = 0;
        goto cleanup;
    }

    /* Handle n even. */
    if (mpz_divisible_ui_p(n, 2)) {
        if (mpz_cmp_ui(n, 2) == 0)
            result = 1;
        else
            result = 0;
        goto cleanup;
    }

    /* Check GCD */
    mpz_mul(res, zD, q);
    mpz_mul_ui(res, res, 2);
    mpz_gcd(res, res, n);
    if ((mpz_cmp(res, n) != 0) && (mpz_cmp_ui(res, 1) > 0)) {
        result = 0;
        goto cleanup;
    }

    /* index = n-(D/n), where (D/n) is the Jacobi symbol */
    int ret;
    mpz_set(index, n);
    ret = mpz_jacobi(zD, n);
    if (ret == -1)
        mpz_add_ui(index, index, 1);
    else if (ret == 1)
        mpz_sub_ui(index, index, 1);

    /* mpz_lucasumod(res, p, q, index, n); */
    mpz_set_si(uh, 1);
    mpz_set_si(vl, 2);
    mpz_set   (vh, p);
    mpz_set_si(ql, 1);
    mpz_set_si(qh, 1);
    mpz_set_si(tmp,0);
    mp_bitcnt_t s;
    s = mpz_scan1(index, 0);
    size_t j;
    for (j = mpz_sizeinbase(index,2)-1; j >= s+1; j--) {
        /* ql = ql*qh (mod n) */
        mpz_mul(ql, ql, qh);
        mpz_mod(ql, ql, n);
        if (mpz_tstbit(index,j) == 1) {
            /* qh = ql*q */
            mpz_mul(qh, ql, q);

            /* uh = uh*vh (mod n) */
            mpz_mul(uh, uh, vh);
            mpz_mod(uh, uh, n);

            /* vl = vh*vl - p*ql (mod n) */
            mpz_mul(vl, vh, vl);
            mpz_mul(tmp, ql, p);
            mpz_sub(vl, vl, tmp);
            mpz_mod(vl, vl, n);

            /* vh = vh*vh - 2*qh (mod n) */
            mpz_mul(vh, vh, vh);
            mpz_mul_si(tmp, qh, 2);
            mpz_sub(vh, vh, tmp);
            mpz_mod(vh, vh, n);
        }
        else {
            /* qh = ql */
            mpz_set(qh, ql);

            /* uh = uh*vl - ql (mod n) */
            mpz_mul(uh, uh, vl);
            mpz_sub(uh, uh, ql);
            mpz_mod(uh, uh, n);

            /* vh = vh*vl - p*ql (mod n) */
            mpz_mul(vh, vh, vl);
            mpz_mul(tmp, ql, p);
            mpz_sub(vh, vh, tmp);
            mpz_mod(vh, vh, n);

            /* vl = vl*vl - 2*ql (mod n) */
            mpz_mul(vl, vl, vl);
            mpz_mul_si(tmp, ql, 2);
            mpz_sub(vl, vl, tmp);
            mpz_mod(vl, vl, n);
        }
    }
    /* ql = ql*qh */
    mpz_mul(ql, ql, qh);

    /* qh = ql*q */
    mpz_mul(qh, ql, q);

    /* uh = uh*vl - ql */
    mpz_mul(uh, uh, vl);
    mpz_sub(uh, uh, ql);

    /* vl = vh*vl - p*ql */
    mpz_mul(vl, vh, vl);
    mpz_mul(tmp, ql, p);
    mpz_sub(vl, vl, tmp);

    /* ql = ql*qh */
    mpz_mul(ql, ql, qh);

    for (j = 1; j <= s; j++) {
        /* uh = uh*vl (mod n) */
        mpz_mul(uh, uh, vl);
        mpz_mod(uh, uh, n);

        /* vl = vl*vl - 2*ql (mod n) */
        mpz_mul(vl, vl, vl);
        mpz_mul_si(tmp, ql, 2);
        mpz_sub(vl, vl, tmp);
        mpz_mod(vl, vl, n);

        /* ql = ql*ql (mod n) */
        mpz_mul(ql, ql, ql);
        mpz_mod(ql, ql, n);
    }

    /* uh contains our return value */
    mpz_mod(res, uh, n);
    if (mpz_cmp_ui(res, 0) == 0)
        result = 1;
    else
        result = 0;

  cleanup:
    mpz_clear(zD);
    mpz_clear(res);
    mpz_clear(index);
    mpz_clear(uh);
    mpz_clear(vl);
    mpz_clear(vh);
    mpz_clear(ql);
    mpz_clear(qh);
    mpz_clear(tmp);

    return result;
}

int GMPY_mpz_is_selfridge_prp(const mpz_t n)
{
    long d = 5, p = 1, q = 0, max_d = 1000000;
    int jacobi = 0;
    mpz_t zD;
    int result;


    /* Take advantage of the cache of mpz_t objects maintained by GMPY2 to
     * avoid memory allocations. */

    mpz_init(zD);

    if (mpz_cmp_ui(n, 2) < 0) {
        result = 0;
        goto cleanup;
    }

    /* Handle n even. */
    if (mpz_divisible_ui_p(n, 2)) {
        if (mpz_cmp_ui(n, 2) == 0)
            result = 1;
        else
            result = 0;
        goto cleanup;
    }


    mpz_set_ui(zD, d);

    while (1) {
        jacobi = mpz_jacobi(zD, n);

        /* if jacobi == 0, d is a factor of n, therefore n is composite... */
        /* if d == n, then either n is either prime or 9... */
        if (jacobi == 0) {
            if ((mpz_cmpabs(zD, n) == 0) && (mpz_cmp_ui(zD, 9) != 0)) {
                result = 1;
                goto cleanup;
            }
            else {
                result = 0;
                goto cleanup;
            }
        }
        if (jacobi == -1)
            break;

        /* if we get to the 5th d, make sure we aren't dealing with a square... */
        if (d == 13) {
            if (mpz_perfect_square_p(n)) {
                result = 0;
                goto cleanup;
            }
        }

        if (d < 0) {
            d *= -1;
            d += 2;
        }
        else {
            d += 2;
            d *= -1;
        }

        /* make sure we don't search forever */
        if (d >= max_d) {
            printf("appropriate value for D cannot be found in is_selfridge_prp()");
            goto cleanup;
        }

        mpz_set_si(zD, d);
    }

    q = (1-d)/4;
    mpz_t p_z, q_z ,n_z;
    mpz_init(p_z);
    mpz_init(q_z);
    mpz_init(n_z);
    mpz_set_si(p_z,p);
    mpz_set_si(q_z,q);

    /* Since "O" is used, the refcount for n is incremented so deleting
     * temp will not delete n.
     */
    result = GMPY_mpz_is_lucas_prp(p_z, q_z, n);

    goto return_result;

  cleanup:

  return_result:
    mpz_clear(zD);

    return result;
}