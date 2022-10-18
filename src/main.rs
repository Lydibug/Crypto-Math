#[allow(dead_code)]

mod crypto_math;
mod misc_functions;

fn factorize_all()
{
    use crate::crypto_math::crypto_math;
    let mut x: u128 = 0;
    loop
    {
        if crypto_math::is_prime(x)
        {
            print!("Prime: ");
        }
        else
        {
            print!("Composite: ");
        }
        let f = crypto_math::factorize(x);
        println!("{} = {} * {}",x,f,x/f);
        x+=1;
    }
}

fn benchmark()
{
    use crate::crypto_math::crypto_math;
    use rand::Rng;
    use std::time::Instant;
    use std::time::Duration;

    let lim: u128 = 1 << 50;//(1 << 127) + ((1 << 127) - 1);

    let mut ans: u128;
    let mut rng = rand::thread_rng();
    let mut runs = 0;
    let mut total = Duration::ZERO;
    let mut average;// = Duration::ZERO;

    loop
    {
        let x: u128 = rng.gen_range(0..lim);
        runs += 1;
        print!("{}: {}...",runs,x);
        let now = Instant::now();
        ans = crypto_math::eulers_phi(x);
        let elapsed = now.elapsed();
        print!("{} Done!",ans);
        total += elapsed;
        average = total/runs;
        println!(" Average: {:.5?}", average);
    }
}

fn main() 
{
    /*
     * Takes a number from the command line and prints a factorization or
     * if the number is prime
     */
    // For using commandline arguments
    use std::env;
    use crate::crypto_math::crypto_math;

    // Get commandline args
    let args: Vec<String> = env::args().collect();
    
    // Check for the correct number of args
    if args.len() < 2
    {
        factorize_all();
    }
    else
    {
        match args[1].as_str()
        {
            "totient" =>
            {
                let number: u128 = args[2].parse().unwrap();
                let totient = crypto_math::eulers_phi(number);
                println!("phi({}) = {}", number, totient);
            },
            "invert" =>
            {
                let number: u128 = args[2].parse().unwrap();
                let modulus: u128 = args[3].parse().unwrap();
                println!("{}^-1 % {} = {}", number, modulus, crypto_math::mod_inv(number,modulus));
            },
            "factor" =>
            {
                let product: u128 = args[2].parse().unwrap();
                // Check if the number is prime
                if crypto_math::is_prime(product)
                {
                    println!("Prime: {}", product);
                }
                else
                {
                    // Print the factorization
                    let factor = crypto_math::factorize(product);
                    println!("{} = {} * {}", product, factor, product/factor);
                }
            },
            "shanks" =>
            {
                let base: u128 = args[2].parse().unwrap();
                let number: u128 = args[3].parse().unwrap();
                let modulus: u128 = args[4].parse().unwrap();
                let log = crypto_math::shanks(base, number, modulus);
                println!("{}",log);
            },
            "bench" => 
            {
                benchmark()
            },
            "challenge" =>
            {
                let numbers = [ 203, 63373, 3254102921,14974358581872710701, 14708973168512866951,165441675456672998379137992434160183811, 223992306552901139411223820675421253569 ];
                    //1021159643789587;
                    //118349971320947570681513251459901556469;
                    //224641912456943976273584427524131892551;
                for number in numbers
                {
                    let factorized = crypto_math::prime_factorize(number);
                    println!("{} = {:?}",number, factorized);
                }
            },
            "legendre" =>
            {
                let number : u128 = args[2].parse().unwrap();
                let prime  : u128 = args[3].parse().unwrap();
                let legendre = crypto_math::legendre_symbol(number, prime);
                println!("({}/{}) = {} ",number, prime, legendre);
            },
            "jacobi" =>
            {
                let upper : u128 = args[2].parse().unwrap();
                let lower  : u128 = args[3].parse().unwrap();
                let jacobi = crypto_math::jacobi_symbol(upper, lower);
                println!("({}/{}) = {} ",upper, lower, jacobi);
            },
            "kronecker" =>
            {
                let upper : u128 = args[2].parse().unwrap();
                let lower  : u128 = args[3].parse().unwrap();
                let kronecker = crypto_math::kronecker_symbol(upper, lower);
                println!("({}/{}) = {} ",upper, lower, kronecker);
            },
            "prime" =>
            {
                let number : u128 = args[2].parse().unwrap();
                print!("{} ",number);
                match crypto_math::is_prime(number) {
                    true => println!("is prime!"),
                    false => println!("is composite")
                };
            },
            "primes" =>
            {
                loop
                {
                    let number = crypto_math::rand_prime();
                    println!("{}",number);
                }
            },
            "primitive" => 
            {
                let number : u128 = args[2].parse().unwrap();
                let primitive : u128 = crypto_math::find_primitive_root(number);
                println!("{} is a primitive root of {}", primitive, number);
            }
            "primitives" =>
            {
                let number : u128 = args[2].parse().unwrap();
                let primitives : Vec<u128> = crypto_math::find_primitive_roots(number);
                for primitive in primitives
                {
                    print!("{}\t", primitive);
                }
                println!();
                //println!("{:?}", primitives);
            }
            _ =>
            {
                println!("incorrect usage")
            }
        }
    }
}

