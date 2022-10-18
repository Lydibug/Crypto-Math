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
        if root > 1
        {
            let mut start = 1;
            let mut end = ( number >> 1 ) + 1;
            let mut mid;
            let mut midsqr;
            while start <= end
            {
                mid = (start + end) >> 1;
                midsqr = mid * mid;
                if midsqr == number
                {
                    root = mid;
                    break;
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
            if root * root != number
            {
                root += 1;
            }
        }
        return root;
    }
}

