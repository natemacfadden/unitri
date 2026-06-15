/*
   Computation of the number of primitive lattice triangulations
   of a rectangle m x n (m is small, n is big)

   Author: Stepan Orevkov  http://picard.ups-tlse.fr/~orevkov

   Last update:  Feb 10, 2017

*/

#include<stdio.h>
#include<stdlib.h>

#define n 7
#define m 6

#define n1 (n+1)
#define n2 (n+2)

#define DEBUG

int only_mem = 0;

int n_primes=200;
int prime[200]=
   {29443, 29453, 29473, 29483, 29501, 29527, 29531, 29537, 29567,
    29569, 29573, 29581, 29587, 29599, 29611, 29629, 29633, 29641, 29663,
    29669, 29671, 29683, 29717, 29723, 29741, 29753, 29759, 29761, 29789,
    29803, 29819, 29833, 29837, 29851, 29863, 29867, 29873, 29879, 29881,
    29917, 29921, 29927, 29947, 29959, 29983, 29989, 30011, 30013, 30029,
    30047, 30059, 30071, 30089, 30091, 30097, 30103, 30109, 30113, 30119,
    30133, 30137, 30139, 30161, 30169, 30181, 30187, 30197, 30203, 30211,
    30223, 30241, 30253, 30259, 30269, 30271, 30293, 30307, 30313, 30319,
    30323, 30341, 30347, 30367, 30389, 30391, 30403, 30427, 30431, 30449,
    30467, 30469, 30491, 30493, 30497, 30509, 30517, 30529, 30539, 30553,
    30557, 30559, 30577, 30593, 30631, 30637, 30643, 30649, 30661, 30671,
    30677, 30689, 30697, 30703, 30707, 30713, 30727, 30757, 30763, 30773,
    30781, 30803, 30809, 30817, 30829, 30839, 30841, 30851, 30853, 30859,
    30869, 30871, 30881, 30893, 30911, 30931, 30937, 30941, 30949, 30971,
    30977, 30983, 31013, 31019, 31033, 31039, 31051, 31063, 31069, 31079,
    31081, 31091, 31121, 31123, 31139, 31147, 31151, 31153, 31159, 31177,
    31181, 31183, 31189, 31193, 31219, 31223, 31231, 31237, 31247, 31249,
    31253, 31259, 31267, 31271, 31277, 31307, 31319, 31321, 31327, 31333,
    31337, 31357, 31379, 31387, 31391, 31393, 31397, 31469, 31477, 31481,
    31489, 31511, 31513, 31517, 31531, 31541, 31543, 31547, 31567, 31573,
    31583};

long alloc_total=0, alloc_pointers=0;
long max_alloc_total=0;
long alloc_memH[2*m*n+1];

int ****H[2*m*n+1];
int *memH[2*m*n+1];

int s[100], a[100], ind[100], indb[100], axx[100];
long n2pow_array[100]={1,1,1,1};
long *n2pow = &(n2pow_array[5]);
int b_array[100]={0}; int *b = &(b_array[1]);
int sgn_array[100]={-1}; int *sgn = &(sgn_array[1]);
int right_vertex_b_array[100] = {-1}; int *right_vertex_b = &(right_vertex_b_array[1]);
int aa0[3], aa1[3], aa2[3];
int incr_sb[100], left_vertex[100], right_vertex[100];
int incr_ind[100], incr_indb[100], incr_indbb[100][3];
int loc_S[100], sb[100];

int ip, p,i,j,ii,S, t, bit, skip, lv,rv,alv,arv,xx,lenseg,a0,am,am1,am2,ind0,memH0;
int b_max[100], am_min, am_max, am1_min, am1_max, am2_min, am2_max;

long mem_page_size, allo, hh, mem;


int gcd( int a, int b ){int c;
    if(a<0)a=-a; if(b<0)b=-b; if(b==0)return a;
    while( (c=a%b) > 0 ){a=b; b=c;}
    return b;
}

int aa(int i){
  switch(b[i]){
    case 0: return a[i];
    case 1: if( a[i]<n1 ){
                if( loc_S[i]==2 )return a[i]-1;
                else return n1;
            }
            else return axx[i];
    case 2: return n1;
  }
}

void usage(char *s){
  printf("\
Usage: %s i\n\
 i is a number in the range 0,1,...%d\n\
   for conputations mod p[i]\n",s,n_primes-1);
  exit(0);
}

void main( int argc, char *argv[] ){ 

    p = 10000;                     // default value of p

    if( argc >= 2 ){
        if(1 != sscanf(argv[1],"%d",&i) )usage(argv[0]);
        if( i >= n_primes )usage(argv[0]);
        p = prime[i]; 
    }
    
    printf("(* Triangulations of rectangle %dx%d *)\np[%d] = %d; f[%d] = {\n",(int)m,(int)n,i,p,i);


    n2pow[0] = 1;
    for( i=0; i<m; i++ )n2pow[i+1]=n2pow[i]*n2;
    for( i=1; i<m-2; i++ )incr_ind[i] = n2pow[m-3-i];

    for( S=0; S<=m; S++ ){
      alloc_pointers += (mem=n2pow[m-3]*sizeof(int ***));
      if( ! only_mem ) H[S] = (int ****)malloc( mem );
      for( i=0; i<n2pow[m-3]; i++ ){
        alloc_pointers += ( mem=n2*sizeof(int **) );
        if( ! only_mem ) H[S][i] = (int ***)malloc( mem );
        for( am2=0; am2<n2; am2++ ){
          alloc_pointers += ( mem=n2*sizeof(int *) );
          if( ! only_mem ) H[S][i][am2] = (int **)malloc( mem );
        }
      }
    }
    alloc_total = alloc_pointers;

//   mem_page_size = 1000000;

   aa0[2] = aa1[2] = aa2[2] = n+1;
//                                       S  i  m2 m1              S  i  m2 m1 m
   if( ! only_mem ){ memH[0] = &memH0; H[0][0][0][0] = memH[0]; H[0][0][0][0][0] = 1; }

// The main loop

   for( S=1; S<=2*m*n; S++ ){

#ifdef DEBUG
printf("====================== S = %d =======================\n",S);
   printf("(* Currenr/Maximal memory allocation: %7.3lf / %6.2lfg\n (including %6.2lfg for pointers) *)\n",
     ((double)alloc_total)/1e9, ((double)max_alloc_total)/1e9, ((double)alloc_pointers)/1e9  );
//printf("mem: %ld\n", alloc_total);
#endif

     if( S>m   )if( ! only_mem ) H[S] = H[S-m-1];
     if( S>m+1 ){ 
        if( ! only_mem )free( memH[S-m-1] );
        alloc_total -= alloc_memH[S-m-1];
     }

/*==================== compute the page size =============================*/

     allo = 0;
     for( i=0; i<=m; i++ ){ a[i]=s[i]=ind[i]=0; }

     while(1){     /* loop over all (m-3)-tuples (a[1],...,a[m-3]), a[i] <= n */
       ind0 = ind[m-3];
       for( am2=0; am2<=n && (s[m-2]=s[m-3]+2*am2) <= S; am2++ ){ a[m-2]=am2;
         aa2[1] = (aa2[0] = am2)-1;
         for( am1=0;  am1<=n && (s[m-1]=s[m-2]+2*am1) <= S; am1++ ){ a[m-1]=am1;
           aa1[1] = (aa1[0] = am1)-1;
           am_min = S - s[m-1] - n; if(am_min<0)am_min=0;
           am_max = (S - s[m-1])/2; if(am_max>n)am_max=n;
           if( am_min <= am_max ){
                allo += (am_max - am_min + 1);
           }
         } /* end of loop by a[m-1] */
       } /* end of loop by a[m-2] */
       /* pass to next (m-3)-tuple (a[1],...,a[m-3]), a[i] <= n */
       for( j=m-3; j>0; j-- )if( a[j]<n && s[j]<S-1 )break;
       if( j<=0 )break;
       a[j]++; s[j]+=2; ind[j] += incr_ind[j];
       while( ++j < m-2 ){ a[j]=0; s[j]=s[j-1]; ind[j]=ind[j-1]; }
     } /* end of loop by (m-3)-tuples (a[1],...,a[m-3]), a[i] <= n */

     for( t=2; t<(1<<m); t+=2 ){  /* the i-th bit of t is 1 <==> a[i]==n1 */
       left_vertex[1]=0;
       for( i=2; i<=m; i++ ){
           if( t&(1<<(i-1)) )left_vertex[i] = left_vertex[i-1];
           else              left_vertex[i] = i-1;
       }
       right_vertex[m-1]=m;
       for( i=m-2; i>=0; i-- ){
           if( t&(1<<(i+1)) )right_vertex[i] = right_vertex[i+1];
           else              right_vertex[i] = i+1;
       }
       loc_S[0] = right_vertex[0];
       loc_S[m] = m - left_vertex[m]; s[0]=ind[0]=0;
       for( i=1; i<m; i++ ){
           ind[i] = ind[i-1]; s[i] = s[i-1];
           loc_S[i] = right_vertex[i] - (lv=left_vertex[i]);
           if( (t&(1<<i)) == 0 ){
               a[i] = 0;
               if( 0 < lv && lv < i-1 )if( a[lv]==0 ){
                    a[i] = 1; s[i] += loc_S[i]; ind[i] += incr_ind[i];
               } 
           }
           else{ a[i]=n1; ind[i] += n1*incr_ind[i]; }
       }
       ind0 = ind[m-3];

       while(1){                                   // loop over   a[1], ...., a[m-3]

           if( t&(1<<(m-2)) ) am2_min = am2_max = n1;
           else{
             am2_min = 0; am2_max = (S - s[m-3])/loc_S[m-2];
             if(am2_max>n)am2_max=n;
           }
           for( am2=am2_min; am2 <= am2_max; am2++ ){
             if(am2 <= n)if((lv=left_vertex[m-2]) < m-3)if(lv>0)if(gcd(am2-a[lv],m-2-lv)>1)continue; 
             aa2[0]=a[m-2]=am2;
             if( loc_S[m-2] == 2 )aa2[1] = am2-1;

             if( am2<=n ){if( (s[m-2] = s[m-3] + am2*loc_S[m-2]) > S )continue;}
             else s[m-2] = s[m-3];
             if( t&(1<<(m-1)) )am1_min = am1_max = n1;
             else{
               am1_min = 0; am1_max = (S - s[m-2])/loc_S[m-1];
               if(am1_max>n)am1_max=n;
             }
             for( am1=am1_min; am1 <= am1_max; am1++ ){
               if(am1 <= n)if((lv=left_vertex[m-1]) < m-2)if(lv>0)if(gcd(am1-a[lv],m-1-lv)>1)continue; 
               aa1[0]=a[m-1]=am1;
               if( loc_S[m-1] == 2 )aa1[1] = am1-1;
               if( am1<=n ){if( (s[m-1] = s[m-2] + am1*loc_S[m-1]) > S )continue;}
               else s[m-1] = s[m-2];

               am_min = (S - s[m-1] - n*loc_S[0] + loc_S[m] - 1)/loc_S[m]; if(am_min<0)am_min=0;
               am_max = (S - s[m-1]) / (loc_S[0] + loc_S[m]);     if(am_max>n)am_max=n;

               if( am_min <= am_max ){ allo += (am_max - am_min + 1); }
             } /* end of loop by a[m-1] */
           } /* end of loop by a[m-2] */

           /* pass to next (m-3)-tuple (a[1],...,a[m-3]) compatible with t */

   loop1:  for( j=m-3; j>0; j-- )if( a[j]<n && s[j] <= S-loc_S[j] )break;
           if( j<=0 )break;
           a[j]++; s[j] += loc_S[j]; // if( a[j]==n1 )s[j]=s[j-1];
           ind[j] += incr_ind[j];
           if( 0 < (lv=left_vertex[j]) )if(j-lv > 1)if( gcd( a[j]-a[lv], j-lv ) > 1 ){
               if( a[j]<n && s[j] <= S-loc_S[j] ){
                   a[j]++; s[j] += loc_S[j]; // if( a[j]==n1 )s[j]=s[j-1];
                   ind[j] += incr_ind[j];
               }
               else goto loop1;
           }
           while( ++j < m-2 ){
               s[j] = s[j-1];
               if( a[j]==n1 ){ ind[j] = ind[j-1] + n1*incr_ind[j]; continue;}
               ind[j] = ind[j-1];
               a[j]=0; lv=left_vertex[j];
               if(j-lv > 1 ){
                   if( a[lv]==0 ){ a[j]=1; s[j] += loc_S[j]; ind[j] += incr_ind[j]; }
                   else if(a[lv]>1)if( gcd(a[lv],j-lv)>1 ){ a[j]=1; s[j] += loc_S[j]; ind[j] += incr_ind[j]; }
               }
           }
           ind0 = ind[m-3];
       } /* end of loop by (m-3)-tuples (a[1],...,a[m-3]), corresponding to t */
     } /* end of loop by t */

     mem_page_size = allo;

/*==================== the page size is computed ========================*/

     alloc_total += ( alloc_memH[S] = mem = mem_page_size*sizeof(int) );
     if( alloc_total > max_alloc_total ) max_alloc_total = alloc_total;

     if( only_mem )continue;
     else memH[S] = (int *)malloc( mem );

     allo = 0;

     for( i=0; i<=m; i++ ){ a[i]=s[i]=b_max[i]=ind[i]=0; }

//==============================================================================

     while(1){     /* loop over all (m-3)-tuples (a[1],...,a[m-3]), a[i] <= n */
       ind0 = ind[m-3];

       for( i=2; i<m-3; i++ )b_max[i] = ( a[i]>0 ? (2*a[i] == a[i-1] + a[i+1] + 1 ? 2 : 1) : 0);

#ifdef _word_9_10_
//       allo = w*((allo+w-1)/w); oo=0;
//       H[S][ind0] = memH[S] + allo/w;
#endif

       for( am2=0; am2<=n && (s[m-2]=s[m-3]+2*am2) <= S; am2++ ){ a[m-2]=am2;
         b_max[m-3] = ( a[m-3]>0 ? (2*a[m-3] == a[m-4] + am2 + 1 ? 2 : 1) : 0);
         aa2[1] = (aa2[0] = am2)-1;

         for( am1=0;  am1<=n && (s[m-1]=s[m-2]+2*am1) <= S; am1++ ){ a[m-1]=am1;
           b_max[m-2] = ( am2>0 ? (2*am2 == a[m-3] + am1 + 1 ? 2 : 1) : 0);
           aa1[1] = (aa1[0] = am1)-1;
           am_min = S - s[m-1] - n; if(am_min<0)am_min=0;
           am_max = (S - s[m-1])/2; if(am_max>n)am_max=n;
           if( am_min <= am_max ){
#ifdef _word_9_10_
//                offs[S][ind0][m-2][m-1] = oo - am_min; oo += (am_max - am_min + 1);
#else
                H[S][ind0][am2][am1] = memH[S] + allo - am_min;
#endif
                allo += (am_max - am_min + 1);
           }
           for( am=am_min; am<=am_max; am++ ){
             a0=a[0] = S - s[m-1] - (a[m]=am);
             aa0[1] = (aa0[0] = am)-1;
             b_max[m-1] = ( am1 >0 ? (2*am1  == am2  + am + 1 ? 2 : 1) : 0);
             b_max[  1] = ( a[1]>0 ? (2*a[1] == a[2] + a0 + 1 ? 2 : 1) : 0);
             b_max[0] = (a0>0 ? 1 : 0);
             b_max[m] = (am>0 ? 1 : 0);

             hh = 0;

             for( i=0; i<=m; i++ ){
                 b[i]=sb[i]=indb[i]=0; sgn[i]=-1;
                 incr_indb[i] = (n+2 - a[i])*incr_ind[i];
             }

             while(1){      /* loop over all (m+1)-tuples (b[0],...,b[m]): THE MOST REPEATED LOOP !!! */

                            /* the array b encodes a subshape of a which we denote by a'              ^/
                            /*  b[i] = 0  means that  a'[i] = a[i]                                    */
                            /*  b[i] = 1  means that  a'[i] = a[i] - 1                                */
                            /*  b[i] = 2  means that  a'[i] = n+1, i.e. the i-th vertex is not in a'  */

                 for( j=m; j>=0; j-- )if( b[j-1]==0 && b[j]<b_max[j] )break;
//printf("trace1: j=%d\n",j);
                 if(j<0)break;
//printf(" indb = ");for( i=0; i<=m-3; i++ )printf("%d ",indb[i]); puts("");

//puts("trace2");
                 if(j>0){
                     if( b[j]++ ){ indb[j] += incr_indb[j]; sb[j]-- ; }
                     else        { indb[j] -= incr_ind[j];  sb[j] += (j==m ? 1:2); }
                 }
                 else{ b[0]++; indb[0] = 0; sb[j]++; }
                 sgn[j]=-sgn[j-1];
                 while( ++j <= m ){ b[j]=0; sb[j]=sb[j-1]; sgn[j]=sgn[j-1]; indb[j]=indb[j-1]; }

                 if( a0 > am || b[0] <= b[m] ){
                     hh += sgn[m] * H[ S-sb[m] ][ ind0 + indb[m-3] ][ aa2[b[m-2]] ][ aa1[b[m-1]] ][ aa0[b[m]] ];
                 }
                 else{
                     for( ii=0,i=m-1; i>2; i-- )ii = ii*n2 + (b[i]==2 ? n1 : a[i]-b[i]);
                     hh += sgn[m] * H[ S-sb[m] ][ii][b[2]==2 ? n1 : a[2]-b[2]][b[1]==2 ? n1 : a[1]-b[1]][a0-1]; 
                 }
             }
             H[S][ind0][am2][am1][am] = (hh%p+p)%p;

#ifdef DEBUG
printf("S=%d  H[%d",S,a0);
for(i=1; i<=m; i++){
  if(a[i]==n1)printf(" ."); else printf(" %d",a[i]);
}
printf("] = %ld\n",hh);
#endif
           } /* end of loop by a[m] */
         } /* end of loop by a[m-1] */
       } /* end of loop by a[m-2] */


       /* pass to next (m-3)-tuple (a[1],...,a[m-3]), a[i] <= n */
       for( j=m-3; j>0; j-- )if( a[j]<n && s[j]<S-1 )break;
       if( j<=0 )break;
       a[j]++; s[j]+=2; ind[j] += incr_ind[j];
       while( ++j < m-2 ){ a[j]=0; s[j]=s[j-1]; ind[j]=ind[j-1]; }

     } /* end of loop by (m-3)-tuples (a[1],...,a[m-3]), a[i] <= n */

     if( S % (2*m) == 0 ){int k=S/2/m;
         ind0=0; for(i=1; i<=m-3; i++)ind0+=n2pow[m-3-i];
         printf("(* f(%d,%d) = *) %d",m,k,H[S][k*ind0][k][k][k]);
         if(k<n)puts(",");else puts("};");
     }
     if( S==2*n*m )break;

//=========================================================================================

     for( t=2; t<(1<<m); t+=2 ){  /* the i-th bit of t is 1 <==> a[i]==n1 */

       left_vertex[1]=0;
       for( i=2; i<=m; i++ ){
           if( t&(1<<(i-1)) )left_vertex[i] = left_vertex[i-1];
           else              left_vertex[i] = i-1;
       }
       right_vertex[m-1]=m;
       for( i=m-2; i>=0; i-- ){
           if( t&(1<<(i+1)) )right_vertex[i] = right_vertex[i+1];
           else              right_vertex[i] = i+1;
       }
       loc_S[0] = right_vertex[0];
       loc_S[m] = m - left_vertex[m]; s[0]=ind[0]=0;
       for( i=1; i<m; i++ ){
           b_max[i] = 0; ind[i] = ind[i-1]; s[i] = s[i-1];
           loc_S[i] = right_vertex[i] - (lv=left_vertex[i]);
           if( (t&(1<<i)) == 0 ){
               a[i] = 0;
               if( 0 < lv && lv < i-1 )if( a[lv]==0 ){
                    a[i] = 1; s[i] += loc_S[i]; ind[i] += incr_ind[i];
               } 
           }
           else{ a[i]=n1; ind[i] += n1*incr_ind[i]; }
       }
       ind0 = ind[m-3]; aa1[2]=aa2[2]=n1;
       if( loc_S[m]   > 1 )aa0[1]=n1;
       if( loc_S[m-1] > 2 )aa1[1]=n1;
       if( loc_S[m-2] > 2 )aa2[1]=n1;

       while(1){                                   // loop over   a[1], ...., a[m-3]

           if( t&(1<<(m-2)) ) am2_min = am2_max = n1;
           else{
             am2_min = 0; am2_max = (S - s[m-3])/loc_S[m-2];
             if(am2_max>n)am2_max=n;
           }
           for( am2=am2_min; am2 <= am2_max; am2++ ){
             if(am2 <= n)if((lv=left_vertex[m-2]) < m-3)if(lv>0)if(gcd(am2-a[lv],m-2-lv)>1)continue; 
             aa2[0]=a[m-2]=am2;
             if( loc_S[m-2] == 2 )aa2[1] = am2-1;

             if( am2<=n ){if( (s[m-2] = s[m-3] + am2*loc_S[m-2]) > S )continue;}
             else s[m-2] = s[m-3];

             if( t&(1<<(m-1)) )am1_min = am1_max = n1;
             else{
               am1_min = 0; am1_max = (S - s[m-2])/loc_S[m-1];
               if(am1_max>n)am1_max=n;
             }

             for( am1=am1_min; am1 <= am1_max; am1++ ){
               if(am1 <= n)if((lv=left_vertex[m-1]) < m-2)if(lv>0)if(gcd(am1-a[lv],m-1-lv)>1)continue; 
               aa1[0]=a[m-1]=am1;
               if( loc_S[m-1] == 2 )aa1[1] = am1-1;

               if( am1<=n ){if( (s[m-1] = s[m-2] + am1*loc_S[m-1]) > S )continue;}
               else s[m-1] = s[m-2];

//   a[m]*loc_S[m] + a[0]*loc_S[0] = S - s[m-1]
//   a[0] >= a[m] >= 0
//   a[0] = (S - s[m-1] - a[m]*loc_S[m])/loc_S[0]
//   a[m] <= (S - s[m-1] - a[m]*loc_S[m])/loc_S[0] <= n
//   loc_S[0]*a[m] <= S - s[m-1] - a[m]*loc_S[m]  &&  S - s[m-1] - a[m]*loc_S[m] <= n*loc_S[0]
//   (loc_S[0] + loc_S[m])*a[m] <= S - s[m]       &&  S - s[m-1] - n*loc_S[0] <= a[m]*loc_S[m]

               am_min = (S - s[m-1] - n*loc_S[0] + loc_S[m] - 1)/loc_S[m]; if(am_min<0)am_min=0;
               am_max = (S - s[m-1]) / (loc_S[0] + loc_S[m]);     if(am_max>n)am_max=n;

               if( am_min <= am_max ){
#ifdef _word_9_10_
//                offs[S][ind0][m-2][m-1] = oo - am_min; oo += (am_max - am_min + 1);
#else
                  H[S][ind0][am2][am1] = memH[S] + allo - am_min;
#endif
                  allo += (am_max - am_min + 1);
               }
               for( am=am_min; am<=am_max; am++ ){
                 if((lv=left_vertex[m]) < m-1)if(lv>0)if(gcd(am-a[lv],m-lv)>1)continue; 

                 if( (xx = S - s[m-1] - (a[m]=am)*loc_S[m]) % loc_S[0] )continue;
                 if( xx < 0 )continue;
                 a0 = a[0] = xx/loc_S[0];
                 if((rv=right_vertex[0])>1)if(gcd(a[rv]-a0,rv)>1)continue;

                 aa0[1] = (aa0[0] = am)-1;

                 hh = b[0] = b[m] = b_max[0] = b_max[m] = sb[0] = sb[m] = indb[0] = 0;
                 sgn[0] = sgn[m] = -1;
                 right_vertex_b[0] = -1;
                 if( right_vertex[0] == 1   ){ b_max[0] = (a0==0 ? 0 : 1); incr_sb[0] = 1; }
                 if(  left_vertex[m] == m-1 ){ b_max[m] = (am==0 ? 0 : 1); incr_sb[m] = 1; }

                 for( i=1; i<m; i++ ){
                     sgn[i] = -1; right_vertex_b[i] = -1;
                     lenseg = loc_S[i];
                     b[i] = b_max[i] = sb[i] = indb[i] = 0;
                     if( (t&(1<<i)) == 0 ){  // the verex is present at i
                         if( lenseg==2 ){    // there are no long adjacent segments
                             if( a[i] > 0 ){
                                 b_max[i] = 1;
                                 if(i<=m-3) incr_indbb[i][0] = - incr_ind[i];
                                 incr_sb[i] = 2;
                                 if( 2*a[i]-1 == a[i-1] + a[i+1] ){
                                     b_max[i] = 2;
                                     if(i<=m-3) incr_indbb[i][1] = (n+1 - a[i])*n2pow[ m-3-i ];
                                 }
                             }
                         }
                         else{                 // there is a long adj segment
                             alv = a[ lv =  left_vertex[i] ];
                             arv = a[ rv = right_vertex[i] ];
                          /* 
                              det( lenseg             i  -   lv  )
                                 ( a[rv] - a[lv]    a[i] - a[lv] )
                                                                   */                  
                             if( lenseg*(a[i] - alv) - (i-lv)*(arv-alv) == 1 ){
                                 b_max[i] = 1;
                                 incr_indbb[i][0] = (n+1 - a[i])*n2pow[ m-3-i ];
                                 incr_sb[i] = 1;
                             }
                         }
                     }
                     else{                 // no vertex at i
                         alv = a[ lv =  left_vertex[i] ];
                         arv = a[ rv = right_vertex[i] ];
                         xx = (alv*(rv-i) + arv*(i-lv))/lenseg;
                      /* 
                          det( i - lv     lenseg      )
                             ( xx-a[lv]   a[rv]-a[lv] )
                                                         */ 
                         if( (i-lv)*(arv-alv) - (xx-alv)*lenseg == 1 ){
                             b_max[i] = 1; axx[i] = xx;
                             if(i<=m-3)incr_indbb[i][0] = (xx - n-1)*n2pow[ m-3-i ];
                             else if(i==m-2)aa2[1] = xx;
                             else if(i==m-1)aa1[1] = xx;
                             incr_sb[i] = 1;
                         }
                     }
                 }

                 while(1){ /* loop over all (m+1)-tuples (b[0],...,b[m]): */

                           /* the array b encodes a subshape of a which we denote by a'             */
                           /*  b[i] = 0  means that  a'[i] = a[i]                                   */
                           /*  b[i] = 1 & a[i] < n1 & no long adjacent segment:   a'[i] = a[i] - 1  */
                           /*  b[i] = 1 & a[i] < n1 &    long adjacent segment:   a'[i] = n1        */
                           /*  b[i] = 1 & a[i] = n1 :                             a'[i] < n1        */
                           /*  b[i] = 2 (appears only if a[i]<n1 & no adj seg):   a'[i] = n1        */

                     for( j=m; j>=0; j-- )if( right_vertex_b[j-1]<j && b[j]<b_max[j] )break;
                     if(j<0)break;
                     right_vertex_b[j] = right_vertex[j];
                     if( j<=m-3 )if( j>0 ) indb[j] = indb[j-1] + incr_indbb[j][ b[j] ];
                     sb[j] = sb[j-1] + (b[j]++ ? 1 : incr_sb[j] );
                     sgn[j]=-sgn[j-1];
                     while( ++j <= m ){
                         right_vertex_b[j] = right_vertex_b[j-1]; 
                         b[j]=0; sb[j]=sb[j-1]; sgn[j]=sgn[j-1]; indb[j]=indb[j-1];
                     }
                     if( a0 > am || b[0] <= b[m] ){
                         hh += sgn[m] * H[ S-sb[m] ][ ind0+indb[m-3] ][ aa2[b[m-2]] ][ aa1[b[m-1]] ][ aa0[b[m]] ];
                     }
                     else{
                         for( ii=0,i=m-1; i>2; i-- )ii = ii*n2 + aa(i);
                         hh += sgn[m] * H[S-sb[m]] [ii] [aa(2)] [aa(1)] [a0-1];
                     }
                 }     // end of loop over b

                 H[S][ind0][am2][am1][am] = (hh%p+p)%p;

#ifdef DEBUG
printf("S=%d  H[%d",S,a0);
for(i=1; i<=m; i++){
  if(a[i]==n1)printf(" ."); else printf(" %d",a[i]);
}
printf("] = %ld\n",hh);
#endif
               } /* end of loop by a[m] */
             } /* end of loop by a[m-1] */
           } /* end of loop by a[m-2] */

           /* pass to next (m-3)-tuple (a[1],...,a[m-3]) compatible with t */

    loop:  for( j=m-3; j>0; j-- )if( a[j]<n && s[j] <= S-loc_S[j] )break;
           if( j<=0 )break;
           a[j]++; b_max[j]=1; s[j] += loc_S[j]; // if( a[j]==n1 )s[j]=s[j-1];
           ind[j] += incr_ind[j];
           if( 0 < (lv=left_vertex[j]) )if(j-lv > 1)if( gcd( a[j]-a[lv], j-lv ) > 1 ){
               if( a[j]<n && s[j] <= S-loc_S[j] ){
                   a[j]++; b_max[j]=1; s[j] += loc_S[j]; // if( a[j]==n1 )s[j]=s[j-1];
                   ind[j] += incr_ind[j];
               }
               else goto loop;
           }
           while( ++j < m-2 ){
               s[j] = s[j-1];
               if( a[j]==n1 ){ ind[j] = ind[j-1] + n1*incr_ind[j]; continue;}
               ind[j] = ind[j-1];
               a[j]=b_max[j]=0; lv=left_vertex[j];
               if(j-lv > 1 ){
                   if( a[lv]==0 ){ a[j]=1; s[j] += loc_S[j]; ind[j] += incr_ind[j]; }
                   else if(a[lv]>1)if( gcd(a[lv],j-lv)>1 ){ a[j]=1; s[j] += loc_S[j]; ind[j] += incr_ind[j]; }
               }
           }
           ind0 = ind[m-3];
       }
     }
   } /* end of loop by S */

   printf("(* Maximal memory allocation: %6.2lfg\n (including %6.2lfg for pointers) *)\n",
     ((double)max_alloc_total)/1e9, ((double)alloc_pointers)/1e9  );
}
