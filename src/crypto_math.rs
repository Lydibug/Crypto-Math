#[allow(dead_code)]
pub mod crypto_math
{

    /*
     * Generates a uniformly random 128 bit number between min and max
     */
    pub fn rand_num
    (min : u128, max : u128) -> u128
    {
        extern crate rand;
        use rand_core::{RngCore, OsRng};
        let mut bytes = [ 0u8; 16 ];
        
        // Fill the bytes with data from the operating systems urandom source
        OsRng.fill_bytes(&mut bytes);
        
        // Convert the bytes to a u128 number and make sure it's less than max
        let number = u128::from_be_bytes(bytes) % max;
        
        // Make sure the number is greater than min
        if number < min
        {
            return min + number;
        }
    
        return number;
    }

    /*
     * Binary GCD algorithm
     */
    pub fn gcd
    (mut u: u128, mut v: u128) -> u128 
    {
        use std::cmp::min;
        use std::mem::swap;
        if u == 0 
        {
            return v;
        } 
        else if v == 0 
        {
            return u;
        }
        let i = u.trailing_zeros();  u >>= i;
        let j = v.trailing_zeros();  v >>= j;
        let k = min(i, j);

        loop {
            debug_assert!(u & 1 == 1, "u = {} is even", u);
            debug_assert!(v & 1 == 1, "v = {} is even", v);

            if u > v 
            {
                swap(&mut u, &mut v);
            }
            v -= u;
            if v == 0 
            {
                return u << k;
            }
            v >>= v.trailing_zeros();
        }
    }

    /*
     * extended euclidian algorithm
     * This can have negative values for s and t
     */
    pub fn gcd_ext
    (a : i128, b : i128) -> (i128, i128, i128)
    {
        let mut a0 : i128 = a;
        let mut b0 : i128 = b;
        let mut t0 : i128 = 0;
        let mut t  : i128 = 1;
        let mut s0 : i128 = 1;
        let mut s  : i128 = 0;
        let mut q  : i128 = a0 / b0;
        let mut r : i128  = a0 - ( q * b0 );
        while r > 0
        {
            let mut temp = t0 - ( q * t );
            t0 = t;
            t = temp;
            temp = s0 - (q * t);
            s0 = s;
            s = temp;
            a0 = b0;
            b0 = r;
            q = a0 / b0;
            r = a0 - (q * b0);
        }
        r = b0;
        return (r, s, t);
    }

    /*
     * Modular multiplication
     * ADD MONTGOMERY REDUCTION
     */
    pub fn mod_mul
    (a: u128, mut b: u128, n: u128) -> u128
    {
        let mut x = 0;
        let mut y = a % n;
        while b > 0
        {
            if b & 1 == 1
            {
                // This line will cause it to fail with numbers >= 2**127
                x = ( x  + y ) % n;
            }
            y = (y << 1) % n;
            b >>= 1;
        }
        return x % n;
    }

    /*
     *  Modular exponentiation using square multiply method
     */
    pub fn mod_pow 
    (mut base: u128, mut exponent: u128, modulus: u128) -> u128
    {
        let mut result: u128 = 1;
        // x**y % 1 == 0 for any x or y
        if modulus <= 1
        {
            return 0;
        }
        
        // Reduce base to a congruent element in the given modulus
        base = base % modulus;
        // If the base is congruent to 0 or 1 in the given modulus, the 
        // result of exponentiation will always be the number itself
        if base <= 1
        {
            return base;
        }

        // Square multiply algorithm
        while exponent > 0
        {
            if exponent & 1 == 1
            {
                result = mod_mul(result,base,modulus);
            }
            exponent >>= 1;
            base = mod_mul(base,base,modulus);
        }
        return result;
    }

    /*
     * Miller-Rabbin primality test
     */
    pub fn is_prime
    (number: u128) -> bool
    {
        extern crate rand;
        use rand::Rng;

        // Use the fastest method for whatever size the number is
        match number
        {
            // Use a lookup table for small numbers
            0 | 1 | 4 | 6 | 8 | 9 | 10 => return false,
            2 | 3 | 5 | 7 | 11 => return true,
            // for numbers 12 to 64, use a custom boolean expression
            12..=64 =>
            {
                if number & 1 == 0
                {
                    return false;
                }
                let bits: [bool; 5] = [ number & 32 == 32, 
                                        number & 16 == 16,
                                        number & 8 == 8,
                                        number & 4 == 4,
                                        number & 2 == 2 ];
                // Boolean expression for some 6 bit even numbers and all 6 bit primes above 10
                // Numbers which wouldn't reach this point were considered dontcares to simplify
                // the expression
                return (!bits[0] && !bits[2] && !bits[3]) ||
                    (!bits[0] && !bits[1] && !bits[4]) ||
                    (!bits[1] &&  bits[2] && !bits[3]) ||
                    (!bits[0] &&  bits[1] &&  bits[3] &&  bits[4]) ||
                    ( bits[1] &&  bits[2] &&  bits[3] && !bits[4]) || 
                    ( bits[0] && !bits[2] &&  bits[3] && !bits[4]) || 
                    ( bits[0] && !bits[1] &&  bits[2] &&  bits[4]) || 
                    ( bits[0] &&  bits[2] && !bits[3] &&  bits[4]); 
            },
            // For numbers greater than 64 bits use Miller-Rabbin
            _ => 
            {
                // Make sure the number ends in 1, 3, 7 or 9
                if number & 1 == 0 || number % 5 == 0
                {
                    return false;
                }

                // Miller-Rabbin primality test with 40 rounds
                // This can be replaced with AKS
                let mut rng = rand::thread_rng();
                let r = (number - 1).trailing_zeros();
                let d =  (number - 1) >> r;

                'WitnessLoop: for _k in 0..40
                {
                    let a = rng.gen_range(2..number - 2);
                    let mut x = mod_pow(a,d,number);
                    if x == 1 || x == number - 1
                    {
                        continue 'WitnessLoop;
                    }
                    for _i in 0.. r - 1
                    {
                        x = mod_mul(x,x,number);
                        if x == number - 1
                        {
                            continue 'WitnessLoop;
                        }
                    }
                    return false;
                }
                return true;
            }
        }
    }

    /*
     * Generates a random bits-size prime number
     * this function is cryptographically secure.
     * This does not ensure the primes will be
     * suitable for use in RSA
     */
    pub fn rand_prime
    (bits : usize) -> u128
    {
        let max = 1 << (bits - 1);
        let min = 2;
        let mut number = rand_num(min, max);
        while !is_prime(number)
        {
            number = rand_num(min, max);
        }
        return number;
    }

    /*
     *  Pollard-Brent factorization
     *  Failure returns None
     */
    pub fn pollard_brent
    (n: u128) -> Option<u128>
    {
        use std::cmp::min;
        // Declare variables
        let mut x : u128;
        let mut y : u128;
        let mut k : u128;
        let mut q : u128;
        let mut d : u128;
        let mut r : u128;
        let mut last_y: u128;
        let m : u128;
        let c : u128;

        y = rand_num(1,n);
        c = rand_num(1,n);
        m = rand_num(1,n);

        last_y = y;
        x = y;
        q = 1;
        r = 1; 
        d = 1;

        // Try pollard brent
        while d == 1
        {
            x = y;
            for _ in 0..r
            {
                y = (mod_mul(y, y, n) + c) % n;
            }
            k = 0;
            while k < r && d == 1
            {
                last_y = y;
                for _ in 0..min(m,r-k)
                {
                    y = (mod_mul(y, y, n) + c) % n;
                    q = mod_mul( if x > y { x - y } else { y - x }, q, n);
                }
                d = gcd(q, n);
                k += m;
            }
            r <<= 1;
        }

        // if it fails, default to pollard rho
        if d == n
        {
            'Backup: loop
            {
                last_y = (mod_mul(last_y, last_y, n) + c) % n;
                d = gcd(if x > last_y { x - last_y } else { last_y - x }, n);
                if d > 1
                {
                    break 'Backup;
                }
            }
        }
        return if d == n { None } else { Some(d) };
    }

    /*
     * Pollard p-1 factorization algorithm
     * Failure returns None
     */
    pub fn pollard_p1
    (number : u128, bound : u128) -> Option<u128>
    {
        let mut a : u128 = 2;
        for j in 2..bound
        {
            a = mod_pow(a, j, number);
        }
        let d =  if a != 0 {  gcd(a - 1, number) } else { 1 };
        
        return if d > 1 && d < number { Some(d) } else { None };
    }

    /*
     * Returns a factor of a given number
     * Failure returns the original number
     */
    pub unsafe fn factorize
    (number: u128) -> u128
    {
        use std::thread;
        use std::sync::mpsc;
        use stop_thread::kill_thread_graceful;
        let factor : u128;
        
        // Attempt pollard p-1 factorization with 
        
        match pollard_p1(number, 1000)
        {
            None => {
                // add a way to close threads when the first thread finds a result
                let (tx, rx) = mpsc::channel();
                let mut handles = vec![];

                for _thread_count in 0..14
                {
                    let tx_clone = tx.clone();
                    let process = thread::spawn( move || 
                        {
                            let result = pollard_brent(number);
                            match tx_clone.send(result)
                            {
                                Ok(_) | Err(_) => { drop( tx_clone ); }
                            };
                        });
                    handles.push(process);
                }
                let result = rx.recv();
                match result
                {
                    Ok(val) => { 
                        factor = match val
                        {
                            None => number,
                            Some(num) => num
                        };
                        drop(tx);
                        for handle in handles
                        {
                            kill_thread_graceful(handle);
                        }
                    },
                    Err(_err) => { factor = number; }
                };
            },
            Some(number) => {
                factor = number;
            }
        }
        /* 
        factor =  match pollard_p1(number, 1000)
            {
                // If pollard p-1 failed, try pollard brent
                // TODO: Add multithreading and a Quadratic Seive
                None => match pollard_brent(number) 
                        { 
                            Some(num) => num,
                            //return the original number if no factor was found
                            None => number 
                        },
                Some(num) => num
            };
        */
        return factor;
    }

    /*
     * Returns a list of the prime factors of a number and their exponants
     * p**k -> (p,k)
     */
    pub unsafe fn prime_factorize
    (mut number : u128 ) -> Vec<(u128,u32)>
    {
        let mut prime_factors : Vec<(u128,u32)> = Vec::<(u128,u32)>::new();
        while number > 1
        {
            //Replace the factorization alg with a quadratic seive for integers <30 digits
            let mut exp = 0;
            let mut factor = factorize(number);
            // Find a prime factor
            while !is_prime(factor)
            {
                factor = factorize(factor);
            }
            // Count how many times that factor appears
            while number % factor == 0
            {
                exp += 1;
                number /= factor;
            }
            prime_factors.push((factor,exp));
        }
        return prime_factors
    }

    /*
     * Eulers Phi function
     */
    pub unsafe fn eulers_phi
    (mut number: u128) -> u128
    {
        let mut pot : u128 = 1;
        let mut phi : u128 = 1;

        if number == 0
        {
            return 1;
        }

        // increases efficiency for even numbers
        let zeros = number.trailing_zeros();
        if  zeros > 0
        {
            let pk : u128 = 1 << zeros;
            let pm = pk >> 1;
            number >>= zeros;
            pot = pk - pm;
        }

        if number == 1
        {
            return pot;
        }

        // increases efficiency for prime numbers
        if is_prime ( number ) 
        {
            phi = number - 1;
        }
        else
        {
            // Make a list of all prime factors of the number
            let mut prime_factors : Vec<(u128,u32)> = prime_factorize(number);
            
            // Loop through each prime factor, using phi(p^k) = p^(k-1) * (p - 1)
            while !prime_factors.is_empty()
            {
                let factor = prime_factors.pop().unwrap();
                let prime = factor.0;
                let exp = factor.1;
                phi *= prime - 1;
                if exp > 1
                {
                    phi *= prime.pow(exp - 1);
                }
            }
        }
        return pot * phi;
    }

    /*
     * Modular square root
     */
    pub fn mod_sqrt
    (number : u128, modulus : u128) -> Option<u128>
    {
        let congruance : u128 = modulus % 4;
        let exponant : u128;
        let root : Option<u128>;
        if congruance == 3
        {
            exponant = (modulus + 1) / 4;
            root = Some(mod_pow(number, exponant, modulus));
        }
        else
        {
            root = None;
        }
        return root;
    }

    /*
     * Determines if a n element of a modulus 
     * is a quadratic residue
     */
    pub fn is_quadratic_residue
    (number : u128, modulus : u128) -> bool
    {
        let exponant : u128 = (modulus - 1) >> 1;
        let residue_check = mod_pow(number, exponant, modulus);
        return residue_check == 1;
    }

    /*
     * Finds the Legendre symbol for two numbers
     * this wont fail or return any warning 
     * if a non-prime or even number is given for p,
     * but will not give a correct result
     */
    pub fn legendre_symbol
    (a : u128, p : u128) -> i8
    {
        if a % p == 0
        {
            return 0;
        }
        else if is_quadratic_residue(a,p)
        {
            return 1;
        }
        else
        {
            return -1;
        }
    }

    /*
     * Finds the jacobi symbol for two numbers
     * this wont fail or return any warning if 
     * an even value of n is given but will not 
     * give a correct result
     */
    pub unsafe fn jacobi_symbol
    (a : u128, n : u128) -> i8
    {
        let n_factors = prime_factorize(n);
        let mut jacobi = 1;

        for p in n_factors
        {
            let prime = p.0;
            let exp = p.1;
            let mut legendre = legendre_symbol(a, prime);
            if exp & 1 == 0
            {
                legendre *= legendre;
            }
            jacobi *= legendre;
        }
        return jacobi;
    }

    /*
     * Calculates the kronecker symbol for two numbers
     */
    pub unsafe fn kronecker_symbol
    (a : u128, n : u128) -> i8
    {
        let n_factors = prime_factorize(n);
        let mut kronecker = 1;

        for p in n_factors
        {
            let prime = p.0;
            let exp = p.1;

            let mut legendre = legendre_symbol(a,prime); 

            if prime == 2
            {
                if a & 1 == 0
                {
                    legendre = 0;
                }
                else
                {
                    let temp = a % 8;
                    if temp == 1 || temp == 7
                    {
                        legendre = 1;
                    }
                    else if temp == 3 || temp == 5
                    {
                        legendre = -1;
                    }
                }
            }

            if exp & 1 == 0
            {
                legendre *= legendre;
            }
            
            kronecker *= legendre;
        }
        return kronecker;
    }

    /*
     * Calculates the inverse of a member of
     * a integer field
     */
    pub unsafe fn mod_inv
    (number: u128, modulus: u128) -> u128
    {
        if number == 1 || number == 0
        {
            return number;
        }
        let phi = eulers_phi(modulus);
        let inv_exp = phi - 1;
        let inverse = mod_pow ( number, inv_exp, modulus );
        return inverse;
    }

    /*
     * Finds a random primitive of a prime field
     */
    pub unsafe fn find_primitive_root
    (number : u128) -> u128
    {
        let totient : u128 = eulers_phi(number);
        let prime_factors : Vec<(u128,u32)> = prime_factorize(totient);
        let mut powers : Vec<u128> = Vec::<u128>::new();
        
        for factor in &prime_factors
        {
            powers.push(totient/ factor.0);
        }
        
        // pick a random elements and check if they're primitive until one is found
        loop
        {
            let primitive = rand_num(2,totient);
            if !is_quadratic_residue(primitive, number)
            {
                let mut is_primitive = true;
                for power in &powers
                {
                    let test = mod_pow(primitive, *power, number);
                    if test == 1
                    {
                        is_primitive = false;
                    }
                }
                if is_primitive
                {
                    return primitive;
                }
            }
        }
    }

    /*
     * Finds all primitive roots of a prime field
     * Returns a vector of primitives in order from least to greatest
     */
    pub unsafe fn find_primitive_roots
    (number : u128) -> Vec<u128>
    {
        let totient : u128 = eulers_phi(number);
        let root_count : u128 = eulers_phi(number - 1);
        let prime_factors : Vec<(u128,u32)> = prime_factorize(totient);
        let mut powers : Vec<u128> = Vec::<u128>::new();
        let mut primitives : Vec<u128> = Vec::<u128>::new();
        
        for factor in &prime_factors
        {
            powers.push(totient/ factor.0);
        }
        
        // itterate over all elements and check if they're primitive until they are all found
        for primitive in 1..number
        {
            // A quadratic residue cannot be a primitive element
            if !is_quadratic_residue(primitive, number)
            {
                let mut is_primitive = true;
                for power in &powers
                {
                    let test = mod_pow(primitive, *power, number);
                    if test == 1
                    {
                        is_primitive = false;
                    }
                }
                if is_primitive
                {
                    primitives.push(primitive);
                }
                if root_count == primitives.len().try_into().unwrap()
                {
                    break;
                }
            }
        }
        return primitives;
    }

    /*
     * Shanks algoritm for the discrete log problem
     * Returns 0 for non-prime moduli
     */
    pub unsafe fn shanks
    (alpha: u128, beta: u128, modulus: u128) -> u128
    {

        use std::collections::HashMap;
        use crate::misc_functions::misc;
        
        if is_prime(modulus)
        {
            let mut l1: Vec<(u128 , u128)> = Vec::<(u128 , u128)>::new();
            let mut l2 = HashMap::new();
            
            // Calculate the ceiling of the square root
            let root = misc::ceil_sqrt(modulus);
            
            // pre-calculate alpha^-1 for use later
            let alpha_inv = mod_inv(alpha, modulus);
            
            // compute (j, a ^ mj mod n) for l1 and (j, b * a ^ -j mod n) for l2
            for j in 0..root
            {
                let alpha_inv_j = mod_pow(alpha_inv, j, modulus);
                l1.push((j, mod_pow(alpha, root * j, modulus)));
                l2.insert( mod_mul(alpha_inv_j, beta, modulus), j);
            }
            // Find an element of each list which have the same second value
            for pair in l1
            {
                let twin = l2.get(&pair.1);
                if !twin.is_none()
                {
                    return mod_mul(root, pair.0, modulus) + twin.unwrap() % modulus;
                }
            }
        }
        return 0;
    }

}
    
