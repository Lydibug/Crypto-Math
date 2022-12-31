#[allow(dead_code)]

pub mod misc
{
    /*
     * Calculates the hamming weight of a number
     */
    pub fn hamming_weight
    (number : u128) -> u32
    {
        let mut weight = 0;
        let mut bits = number;
        while bits > 0
        {
            weight += bits & 1;
            bits >>= 1;
        }
        return weight.try_into().unwrap();
    }

    /*
     * Calculates the hamming distance between two numbers
     */
    pub fn hamming_dist
    (num_a : u128, num_b : u128) -> u32
    {
        let diff = num_a ^ num_b;
        return hamming_weight(diff);
    }

    /*
     * Returns a number in NAF form
     * returned value is given as a string of positive bits and a string of
     * negative bits (pos, neg)
     */
    pub fn to_naf
    (number : u128) -> (u128, u128)
    {
        let num_half = number >> 1;
        let num3 = number + num_half;
        let c = num_half ^ num3;
        return (num3 & c, num_half & c);
    }

    /*
     * Takes a number in NAF and returns it as a binary integer
     */
    pub fn from_naf
    ( ( np, nm ) : (u128, u128) ) -> u128
    {
        return np - nm;
    }

    /*
     * Calculates the ceiling of the square root of a number
     * Uses a binary search to find the square root
     *
     * Fails on numbers over 64 bits
     */
    pub fn ceil_sqrt
    (number: u128) -> u128
    {

        let mut root = number;
        // if number == 1 or 0, sqrt(number) = number and no further
        // calculation is needed
        if root > 1
        {
            let mut start = 1;
            // Skip the first itteration of the search since 
            // Sqrt(x) !> x/2
            let mut end = ( number >> 1 ) + 1;
            let mut mid;
            let mut midsqr;
            // Begin binary search
            while start <= end
            {
                mid = (start + end) >> 1;
                midsqr = mid * mid;
                // If the number is a perfect square
                // and the square root is found
                if midsqr == number
                {
                    return mid;
                }
                if midsqr < number
                {
                    start = mid + 1;
                    root = mid;
                }
                else
                {
                    end = mid - 1;
                }
            }
            // If the number is not a perfect square,
            // then root < sqrt(number) so add one to the 
            // root
            if root * root != number
            {
                root += 1;
            }
        }
        return root;
    }

    // Square multiply exponentiation
    fn fast_exp
    (a : u128, b : u32) -> u128
    {
        let mut base = a;
        let mut exp = b;
        let mut ans = 1;

        while exp > 0
        {
            if exp & 1 == 1
            {
                ans = ans * base;
            }
            base = base * base;
            exp >>= 1;
        }
        return ans;
    }

    /*
     * Checks for a perfect power and returns the root
     */
    pub fn perfect_power
    (number : u128) -> Option<u128>
    {
        use::std::sync::{Arc,Mutex};
        use::std::thread;

        let root = Arc::new(Mutex::new(0));
        let mut handles = vec![];

        if number == 1
        {
            return Some(1);
        }

        // Handles powers of two
        if hamming_weight(number) == 1
        {
            return Some(2);
        }
        
        // Find log base 2 of n
        let log_2n = 128 - number.leading_zeros() - 1;
        for power in 2 .. log_2n // Check powers concurrently
        {
            let root = Arc::clone(&root);
            let handle = thread::spawn( move || 
                {
                    let mut low_base = 1;
                    let mut high_base = 1 << (log_2n / power + 1 );
                    // Binary search for a value for which base ** power = number
                    while low_base < high_base - 1
                    {
                        let mid_base = (low_base + high_base) >> 1;
                        let base_pow = fast_exp(mid_base, power);
                        if base_pow > number
                        {
                            high_base = mid_base;
                        }
                        else if base_pow < number
                        {
                            low_base = mid_base;
                        }
                        else
                        {
                            let mut rt = root.lock().unwrap();
                            *rt = mid_base;
                            break;
                        }
                    }
                });
            handles.push(handle);
        }

        // Close all threads
        for handle in handles
        {
            handle.join().unwrap();
        }

        let result = *root.lock().unwrap();
        // Failure, return None
        if result > 0
        {
            return Some(result);
        }
        return None;
    }
}

