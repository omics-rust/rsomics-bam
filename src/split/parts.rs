pub(super) struct SplitMix64(u64);

impl SplitMix64 {
    pub(super) fn new(seed: u64) -> Self {
        Self(seed)
    }

    pub(super) fn index(&mut self, upper: usize) -> usize {
        assert!(upper > 0, "SplitMix64 upper bound must be positive");
        let upper = u64::try_from(upper).expect("usize fits in u64 on supported targets");
        usize::try_from(self.next() % upper).expect("bounded SplitMix64 index fits in usize")
    }

    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut value = self.0;
        value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        value ^ (value >> 31)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splitmix64_seed_zero_vector_is_stable() {
        let mut generator = SplitMix64::new(0);
        assert_eq!(
            (0..8).map(|_| generator.next()).collect::<Vec<_>>(),
            [
                0xe220_a839_7b1d_cdaf,
                0x6e78_9e6a_a1b9_65f4,
                0x06c4_5d18_8009_454f,
                0xf88b_b8a8_724c_81ec,
                0x1b39_896a_51a8_749b,
                0x53cb_9f0c_747e_a2ea,
                0x2c82_9abe_1f45_32e1,
                0xc584_133a_c916_ab3c,
            ]
        );
    }
}
