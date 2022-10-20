#[allow(dead_code)]

pub mod misc
{
    /*
     * Calculates the ceiling of the square root of a number
     * Uses a binary search to find the square root
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
}

