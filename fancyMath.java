import java.io.IOException;
/*a*/
class MathFunc {
    
    MathFunc() {
        /*define here initial values for upper, lower and epsilon*/
        this(10000.0, -10000.0, 0.00005);
    }
    
    MathFunc(double upper_t, double lower_t, double epsilon_t) {
        upper = upper_t; 
        lower = lower_t;
        epsilon = epsilon_t;
        runs = 0;
        last_middle = 1.0;
    }
    
    public double upper, lower, epsilon, last_middle;
    public int runs;
    
    static double function_f(double x) {
        return (Math.pow(x, 3) - 24 * Math.pow(x, 2) + 59 * x + 420);
        //return (-1/Math.exp(x) + 1e20);
    }

    double find_zero_rec(double upper_t, double lower_t) {
        
        runs++;     
        
        double middle = lower_t + 0.5 * Math.abs(upper_t - lower_t);
        double middle_y = function_f(middle);

        // another break, there are functions where the zero can't be reached due to insufficient precision of double,
        // probably not sufficient for all kind of functions
        if( Math.abs(last_middle - middle) <= epsilon)
            return middle;      
        last_middle = middle;
        if( middle_y > 0.0 )
             upper_t = middle;
        else if ( middle_y < 0.0)
            lower_t = middle;                 
                    
        if(Math.abs(middle_y - 0.0) <= epsilon)
            return middle;
        else
            return find_zero_rec(upper_t, lower_t);
                 
    }
    
    double find_zero() {
        double zero = Double.NaN;  
        // check for valid bounds
        if( upper != lower && upper > lower) {
            // first, check if we can have zeroes 
            double middle = 0.0;
            double up_res = function_f(upper);
            double lo_res = function_f(lower);
            double up_sign = Math.signum(up_res);
            double lo_sign = Math.signum(lo_res);
        
            // TODO : check with Double.isNaN() for ... 
        
            if( up_sign == 0.0)  // seldom case, if upper or lower already zero
                return upper;
            if( lo_sign == 0.0)
                return lower;
        
            double sign_res = up_sign + lo_sign;
            if( sign_res == 0.0) { // -1.0+1.0 , if different sign and therefore zero possible
            
                boolean recursive = true;
                if(!recursive) {
               
                    double middle_y = 1.0;
                    while( Math.abs(middle_y - 0.0) > epsilon) {
                        middle = lower + 0.5 * Math.abs(upper - lower);
                        runs++;
                        middle_y = function_f(middle);
                        // another break, there are functions where the zero can't be reached due to insufficient precision of double,
                        // probably not sufficient for all kind of functions
                        if( Math.abs(last_middle - middle) <= epsilon)
                            break;      
                        last_middle = middle;
                        if( middle_y > 0.0 )
                            upper = middle;
                        else if ( middle_y < 0.0)
                            lower = middle;                 
                          
                    }
                    zero = middle; 
                } else {
                    zero = find_zero_rec(upper, lower);
                }

            } else if( sign_res == -2.0 || sign_res == 2.0 )  // same sign, no zero possible
                return Double.NaN;
            //else 
            // TODO : handle Inf, NaN .. 

        }
 
        return zero;    
    }
   
}
/*b*/
class MathPower {
    
    MathPower() {
       mults = calls = 0;
    }
    
    public int mults, calls;
    /*only for counting multiplications and calls these methods need to be non-static*/
    /*iterative*/
    double pow_iter(double base, int exponent) {
        double result = 1.0;
        int sign = Integer.signum(exponent);
        exponent = Math.abs(exponent);
        for(int i = 0; i < exponent; i++, mults++)
            result *= base;
        return (sign < 0 ? 1/result : result);        
    }
    
    /*recursive*/
    double pow_recur(double base, int exponent) {
        if(exponent == 0)
            return 1.0;
        if(Integer.signum(exponent) < 0) {
            mults++;
            calls++;
            return (1/(base * pow_recur(base, Math.abs(exponent) - 1))); 
        }
        mults++; 
        calls++;      
        return base * pow_recur(base, exponent - 1);
    }
    /*more efficient recursive*/
    double pow_recur_eff(double base, int exponent) {
        if(exponent == 0)
            return 1.0;
        if(Integer.signum(exponent) < 0) {
            mults += 2;
            calls += 2;
            double tmp = pow_recur(base, Math.abs(exponent)/2 );
            return (1/(pow_recur(base, Math.abs(exponent)%2) * tmp * tmp));  
        }
        mults += 2;  
        calls += 2;
        double tmp2 = pow_recur(base, exponent/2);
        return pow_recur(base, exponent%2) * tmp2 * tmp2;    
    }
}
/*c, probably not a good solution*/
class fancyPrinter2 {
    
    fancyPrinter2() {}

     public static void printPattern(int n){

        if ( n==-1 )
            return;

        printPattern( n-1 );

        for( int i = 0; i <= n; i++ )
            System.out.print( n );

        System.out.println();

        printPattern( n-1 );
    }
   
   
}
/*d + e*/
class PrimeTest {
    
    PrimeTest() {}
    /*optimize: level of optimization 0 == no optimization*/
    /*probably not the best idea for x to be > 2e53 */
    static boolean isPrime(long n, int optimized) {
        
        boolean is_prime = true;

        if ( n <= 1  ) {
            is_prime = false;
        } else if( n == 2 ) {
            is_prime = true;
        } else {
            
            long i = 2, max_div = 0;
            /*f 1)*/
            // optimization: only need to divide till sqrt(n)
            if( optimized == 1 )
                max_div = (long)Math.ceil(Math.sqrt((double)n));
            else if( optimized == 0 )
                max_div = n - 1;
            
            while ( i <= max_div ) {
   /*e 3), -->*/if( ( n % i ) == 0 ) { /*<-- on optimized == 0, this is calculated n - 2 times, on optimized == 1,  ~ 1/2*sqrt(n) times */ 
                    is_prime = false;
                    break;
                }
                /*f 1)*/
                /* another optimization: no need to divide by other even numbers except 2*/
                if( optimized == 1 )
                    i = (i == 2) ? i+1 : i+2; 
                else if( optimized == 0 )
                    i++;
            }           
            
        }    
        return is_prime;
    }
    /*f 2). ??? maybe this is meant, dividing only by primes ?, trade off here of course you have to use hardcoded list,
    there is also upper limit for n, which for this list is 1026169*/
    static int isPrime_list(long n) {
        
        long[] primes_table = {   2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 
                            31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 
                            73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 
                            127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 
                            179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 
                            233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 
                            283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 
                            353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 
                            419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 
                            467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 
                            547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 
                            607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 
                            661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 
                            739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 
                            811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 
                            877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 
                            947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013  };
         
         
        // of course this assumes that your list is sorted 
        if( n > primes_table[primes_table.length - 1] * primes_table[primes_table.length - 1] ) 
            return -1; // n too big
            
        int is_prime = 1;    
   
        if ( n <= 1  ) {
            is_prime = 0;
        } else if( n == 2 ) {
            is_prime = 1;
        } else {
            
            /* I think first searching if n is already is in the primes_table is not worth on average. */
            
          
            long div = (long)Math.ceil(Math.sqrt((double)n)), nmb = 0;

            /*binary search for the biggest divider in the primes table*/
            int upper_bound = primes_table.length, lower_bound = 0, middle = (upper_bound - lower_bound)/2;
            while( upper_bound - lower_bound >= 2 ) {
                nmb = primes_table[middle];
                if( nmb > div)
                    upper_bound = middle;
                else if( nmb < div)
                    lower_bound = middle;
                else
                    break;       
                
                middle = lower_bound + (upper_bound - lower_bound) / 2;
            }
            // primes_table[middle+1] is always bigger than div (equal is catched in the loop above)*/
            if( primes_table[middle] < div )
                middle++;
            
            /*trial division loop*/
            for(int i = 0; i < middle; i++)
                if( n % primes_table[i] == 0 ) {
                    is_prime = 0;
                    break;
                }        
           
            
        }    
        return is_prime;
    }
    
    /*for e 1),    set save to true*/
    /* method == 0 Trial division, method == 1 Trial division with primes Lookup Table, method == 2 Sieve Of Eratosthenes*/
    static long[] findPrimes(long till_nmb, int[] primes_found, boolean save, int method, int optimization) {
        long[] buff = null;
        boolean prime_flag = false;
        if(save)
            buff = new long[(int)till_nmb / 2 + 2]; // we cannot know how much primes we have, but should use normal array ...
        int j = 0;
        for(long i = 2; i < till_nmb;) {
            if( method == 0)
                prime_flag = isPrime(i, optimization);
            else if( method == 1)
                prime_flag = isPrime_list(i) == 1;

            if( prime_flag )
                if(save)
                    buff[j++] = i;
            i++;
            //i = (i == 2) ? i+1 : i+2; another possible optimization here          
        }
        
        primes_found[0] = j;
        return buff;
    }
    /*e 2)*/  
    /* returns array with the single steps, containing times at every step, one step calcs primes till steps * stepSize */
    /* method == 2 Sieve Of Eratosthenes */
    static double[] timeForFindPrimes(int steps, int stepSize, int method) {
        double[] buff = new double[steps];
        long t_beg = 0, t_end = 0;
        long till_nmb = 0;
        int[] primes_nmb = new int[1];
        SieveOfEratosthenes er_sieve = null;
        for(int i = 1; i < steps + 1; i++) {
            till_nmb = i * stepSize;         
            if( method == 2 )
                er_sieve = new SieveOfEratosthenes((int)till_nmb);
            t_beg = System.nanoTime();  
            if( method == 2 )
               er_sieve.sieve();
            else           
                findPrimes(till_nmb, primes_nmb, false, method, 0);
            t_end = System.nanoTime();       
            buff[i - 1] = (t_end - t_beg) / 10e6; // time in milliseconds
        }
        
        return buff;
    }
    /*e 4)*/
    static double[] timeForFindPrimes_exp(int steps, int base, int method) { 
        double[] buff = new double[steps];
        long t_beg = 0, t_end = 0;
        long till_nmb = 0;
        int[] primes_nmb = new int[1];
        SieveOfEratosthenes er_sieve = null;
        for(int i = 0; i < steps; i++) {
            till_nmb = (long)Math.pow(base, i);
            if( method == 2 )
                er_sieve = new SieveOfEratosthenes((int)till_nmb);
            t_beg = System.nanoTime();      
            if( method == 2 )
               er_sieve.sieve();
            else           
                findPrimes(till_nmb, primes_nmb, false, method, 0);
            t_end = System.nanoTime();
            buff[i] = (t_end - t_beg) / 10e6;
        }
        
        return buff;
    }
    
}
/*f 3)*/
class SieveOfEratosthenes {
    
    SieveOfEratosthenes() {
        this(1000000);
    }
    /* create a new field, every element is one number, and initialize all to true, 
    later these elements which remain true are the primes */
    SieveOfEratosthenes(int n) {
        max = n;
        sieve_field = new boolean[max + 1]; // max + 1 for easy conversion of index to number, to be able to begin counting at 1 
        for(int i = 0;i < max + 1; i++)
            sieve_field[i] = true;

    }
    
    boolean[] sieve_field;
    int max;
    
    void sieve() {
        if( max == 0)
            return;
        /* Iterate over every element i (till sqrt(max)) in the sieve_field beginning at i = 2 and if the element i is true, then begin by its first multiple (which is i * 2)
        to mark all its other multiples as false, this way all multiples are eliminated (marked as false) and the remaining elements are the primes. */
        for(int i = 2; i*i <= max; i++)
            if( sieve_field[i] )
                for(int j = i * 2; j <= max; j += i)
                    sieve_field[j] = false;
        
    }
    
    void print_fancy() {
        final int max_per_line = 10;
        int line_cnt = max_per_line;
        for(int i = 2; i <= max; i++)
            if( sieve_field[i] ) 
                if( line_cnt == max_per_line ) {
                    System.out.print(i); 
                    line_cnt--;
                } else if( line_cnt == 0 ) {
                    System.out.print(", "); 
                    System.out.println(); 
                    line_cnt = max_per_line;                     
                } else {
                    System.out.print(", " + i); 
                    line_cnt--; 
                }     
    }          
    
}

public class fancyMath {
    
    public static void printLongArray(long[] array, int elems, int chars_per_line) {
        int cnt = 0;
        int left = elems;
            for(long i = 0; i < elems/chars_per_line; i++) {
                for(long j = 0; j < chars_per_line; j++, left--) 
                    if(left == 1) // last line
                        System.out.print(array[cnt++]);
                    else 
                        System.out.print(array[cnt++] + ", ");
                System.out.println();        
            }
            while( left --> 0 ) {
                if(left == 0) // last line
                    System.out.print(array[cnt++]);
                else 
                    System.out.print(array[cnt++] + ", ");
            }       
    }
    
    public static void printDoubleArray(double[] array, int elems, int chars_per_line) {

        
        int cnt = 0;
        int left = elems;
            for(long i = 0; i < elems/chars_per_line; i++) {
                for(long j = 0; j < chars_per_line; j++, left--) 
                    if(left == 1) // last line
                        System.out.print(array[cnt++]);
                    else 
                        System.out.print(array[cnt++] + ", ");
                System.out.println();        
            }
            while( left --> 0 ) {
                if(left == 0) // last line
                    System.out.print(array[cnt++]);
                else 
                    System.out.print(array[cnt++] + ", ");
            }  
    }
    public static void main(String[] args) {
        /*a*/ boolean A = true;
        if(A) {
        MathFunc ob = new MathFunc();
        double zero = ob.find_zero();
        if( !Double.isNaN(zero) ) {
            System.out.println("Nullstelle gefunden bei x = " + zero + ", Funktionswert f(x) = " + MathFunc.function_f(zero) + " , " + ob.runs +" Überprüfungen durchgeführt.");    
        } else {
            System.out.println("Keine Nulstelle gefunden, Bereich falsch oder keine Nullstelle in dem Bereich vorhanden."); 
        }
        
        try {
            System.in.read();
            System.in.read(); 
        } catch (IOException e) {

        }
        } 
        /*b*/ boolean B = true;
        if(B) {
        MathPower ob2 = new MathPower();
        System.out.println("Potenzierung einer Zahl:");
        System.out.println("Iterativ: 5.0^29 = " + ob2.pow_iter(5.0, 29) + ", 3.0^7 = " + ob2.pow_iter(3.0, 7) + 
        ", 10.0^-2 = " + ob2.pow_iter(10.0, -2) + ", 9.0^0 = " + ob2.pow_iter(9.0, 0)); 
        System.out.println("Iterativ: Multiplikationen ausgeführt: " + ob2.mults);
        ob2.mults = 0; 
        System.out.println("Rekursiv: 5.0^29 = " + ob2.pow_recur(5.0, 29) + ", 3.0^7 = " + ob2.pow_iter(3.0, 7) + 
        ", 10.0^-2 = " + ob2.pow_recur(10.0, -2) + ", 9.0^0 = " + ob2.pow_recur(9.0, 0));
        System.out.println("Rekursiv: Multiplikationen ausgeführt: " + ob2.mults + ", Methodenaufrufe ausgeführt: " + ob2.calls);
        ob2.mults = 0; 
        ob2.calls = 0;        
        System.out.println("Rekursiv, effizienter: 5.0^29 = " + ob2.pow_recur_eff(5.0, 29) + ", 3.0^7 = " + ob2.pow_iter(3.0, 7) + 
        ", 10.0^-2 = " + ob2.pow_recur_eff(10.0, -2) + ", 9.0^0 = " + ob2.pow_recur_eff(9.0, 0));  
        System.out.println("Rekursiv, effizienter: Multiplikationen ausgeführt: " + ob2.mults + ", Methodenaufrufe ausgeführt: " + ob2.calls);
        ob2.mults = 0; 
        ob2.calls = 0; 
        try {
            System.in.read();
            System.in.read();
        } catch (IOException e) {

        }
        }
        
        /*c*/ boolean C = true;
        if(C) {
            System.out.println("Fancy Printer2:");
            fancyPrinter2.printPattern(3); 
            try {
                System.in.read();
                System.in.read();
            } catch (IOException e) {

            }
        }
        
   
        final int until = 15000, stepsize = 1000;
        /*e*/ boolean E = true;
        if(E) {
            final long max = 9960;
            final int chars_per_line = 10;
            /*e 1)*/
            System.out.println("Primzahlen bis " + max + ", unoptimierte Probedivision :");
            int[] primes_cnt = new int[1]; // ...
            long[] buff = PrimeTest.findPrimes( max, primes_cnt, true, 0, 0);
            printLongArray(buff, primes_cnt[0], chars_per_line);
            
            System.out.println();  try {
                System.in.read();
                System.in.read();
            } catch (IOException e) {

            }
            
            /*e 2)*/     
            System.out.println("Laufzeitmessung der unoptimierten Probedivision, bei Obergrenze für Primzahlen von " + until + " mit Schrittweite " + stepsize + ", Zeit in Millisekunden: ");           
            double[] times = PrimeTest.timeForFindPrimes(until/stepsize, stepsize, 0); 
            printDoubleArray(times, times.length, 2); 
            
            /* possible output max 10000, step 1000:
              0.1553151, 
            0.1406001, 
            0.2108051, 
            0.3554401, 
            0.5322582, 
            0.7506107, 
            1.003425, 
            1.3396406, 
            1.5913394, 
            1.9539233
             */
            
            try {
                System.in.read();
                System.in.read();
            } catch (IOException e) {

            }
            
            
            /*e 4)*/            
            System.out.println();
            final int base = 2, exp = 16;
            System.out.println("Laufzeitmessung unoptimierten Probedivision bei exponentiellem Wachstum zur Basis " + base + " mit maximalem Exponent " + exp + ", Zeit in Millisekunden: "); 
            times = PrimeTest.timeForFindPrimes_exp(exp, base, 0);  
            printDoubleArray(times, times.length, 1);   
            
            /* possible output: 
            7.12E-5, 
            2.38E-5, 
            3.085E-4, 
            2.38E-5, 
            2.37E-5, 
            7.12E-5, 
            2.373E-4, 
            6.883E-4, 
            0.0022785, 
            0.0080221, 
            0.0277687, 
            0.098591, 
            0.3775364, 
            1.3163813, 
            4.8656551, 
            17.9855835
             */
             
            try {
                System.in.read();
                System.in.read();
            } catch (IOException e) {

            }
            
            /*f 2)*/
            
            System.out.println();
            System.out.println("Laufzeitmessung Probedivision mit vorberechneter Liste von Primzahlen, bei Obergrenze für Primzahlen von " + until + " mit Schrittweite " + stepsize + ", Zeit in Millisekunden: ");           
            times = PrimeTest.timeForFindPrimes(until/stepsize, stepsize, 1); 
            printDoubleArray(times, times.length, 1);     
            
            /* possible output with max 10000, step 1000:
            0.4279947, 
            0.2462637, 
            0.3797198, 
            0.385867, 
            0.3842294, 
            0.4789279, 
            0.5101143, 
            0.3020385, 
            0.350005, 
            0.4792839
             */
            
            long[] new_primes = PrimeTest.findPrimes(1000000, primes_cnt, true, 1, 0);
            printLongArray(new_primes, primes_cnt[0], chars_per_line);
            
            try {
                System.in.read();
                System.in.read();
            } catch (IOException e) {

            }
            
        }
        
        
        /*Sieve of Eratosthenes*/
        boolean SIEVE_ER = true;
        if( SIEVE_ER ) {
            final int max_er = 1000000;
            SieveOfEratosthenes er_sieve = new SieveOfEratosthenes(max_er);
            er_sieve.sieve();
            System.out.println("Primzahlen bis mit Sieb des Erathostenes " + max_er + " :");
            er_sieve.print_fancy();
            
            try {
                System.in.read();
                System.in.read();
            } catch (IOException e) {

            }
            
            System.out.println();
            System.out.println("Laufzeitmessung beim Sieb des Erathostenes, bei Primzahlen bis " + until + " mit Schrittweite " + stepsize + ", Zeit in Millisekunden: ");           
            double[] times = PrimeTest.timeForFindPrimes(until/stepsize, stepsize, 2); 
            printDoubleArray(times, times.length, 1); 
            /* possible output with 100000 max and 1000 step:
                  0.0079746, 
                8.307E-4, 
                0.0013529, 
                0.0017088, 
                0.0023022, 
                0.0027531, 
                ...
                0.0334174, 
                0.0242324, 
                0.024446, 
                0.0247546
             */
            
            
        }
            
    }
    
}
