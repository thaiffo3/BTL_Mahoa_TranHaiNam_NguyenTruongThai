package ecm;


import java.math.BigInteger;
import java.util.Arrays;

public class SievePrimeNumber {

	private static void enumber1(int f, int g, int d, long L, int B) {
		int x = f, y0 = g;
		long temp = L + B;
		int k0 = ((f << 2) * f + g * g - d) / 60;

		while (k0 < temp) {
			k0 += (x << 1) + 15;
			x += 15;
		}

		while (true) {
			x -= 15;
			k0 -= (x << 1) + 15;
			if (x <= 0)
				break;

			while (k0 < L) {
				k0 += y0 + 15;
				y0 += 30;
			}

			int k = k0, y = y0;
			while (k < temp) {
				segs[d][(int) (k - L) >> 5] ^= 1 << ((k - L) & 31);
				k += y + 15;
				y += 30;
			}
		}
	}

	public static BigInteger[] smallPrimes;

	private static final int[][] DF1 = { { 1, 0, 1 }, { 1, 0, 11 }, { 1, 0, 19 },
			{ 1, 0, 29 }, { 1, 2, 15 }, { 1, 3, 5 }, { 1, 3, 25 }, { 1, 5, 9 },
			{ 1, 5, 21 }, { 1, 7, 15 }, { 1, 8, 15 }, { 1, 10, 9 },
			{ 1, 10, 21 }, { 1, 12, 5 }, { 1, 12, 25 }, { 1, 13, 15 },
			{ 13, 1, 3 }, { 13, 1, 27 }, { 13, 4, 3 }, { 13, 4, 27 },
			{ 13, 6, 7 }, { 13, 6, 13 }, { 13, 6, 17 }, { 13, 6, 23 },
			{ 13, 9, 7 }, { 13, 9, 13 }, { 13, 9, 17 }, { 13, 9, 23 },
			{ 13, 11, 3 }, { 13, 11, 27 }, { 13, 14, 3 }, { 13, 14, 27 },
			{ 17, 2, 1 }, { 17, 2, 11 }, { 17, 2, 19 }, { 17, 2, 29 },
			{ 17, 7, 1 }, { 17, 7, 11 }, { 17, 7, 19 }, { 17, 7, 29 },
			{ 17, 8, 1 }, { 17, 8, 11 }, { 17, 8, 19 }, { 17, 8, 29 },
			{ 17, 13, 1 }, { 17, 13, 11 }, { 17, 13, 19 }, { 17, 13, 29 },
			{ 29, 1, 5 }, { 29, 1, 25 }, { 29, 4, 5 }, { 29, 4, 25 },
			{ 29, 5, 7 }, { 29, 5, 13 }, { 29, 5, 17 }, { 29, 5, 23 },
			{ 29, 10, 7 }, { 29, 10, 13 }, { 29, 10, 17 }, { 29, 10, 23 },
			{ 29, 11, 5 }, { 29, 11, 25 }, { 29, 14, 5 }, { 29, 14, 25 },
			{ 37, 2, 9 }, { 37, 2, 21 }, { 37, 3, 1 }, { 37, 3, 11 },
			{ 37, 3, 19 }, { 37, 3, 29 }, { 37, 7, 9 }, { 37, 7, 21 },
			{ 37, 8, 9 }, { 37, 8, 21 }, { 37, 12, 1 }, { 37, 12, 11 },
			{ 37, 12, 19 }, { 37, 12, 29 }, { 37, 13, 9 }, { 37, 13, 21 },
			{ 41, 2, 5 }, { 41, 2, 25 }, { 41, 5, 1 }, { 41, 5, 11 },
			{ 41, 5, 19 }, { 41, 5, 29 }, { 41, 7, 5 }, { 41, 7, 25 },
			{ 41, 8, 5 }, { 41, 8, 25 }, { 41, 10, 1 }, { 41, 10, 11 },
			{ 41, 10, 19 }, { 41, 10, 29 }, { 41, 13, 5 }, { 41, 13, 25 },
			{ 49, 0, 7 }, { 49, 0, 13 }, { 49, 0, 17 }, { 49, 0, 23 },
			{ 49, 1, 15 }, { 49, 4, 15 }, { 49, 5, 3 }, { 49, 5, 27 },
			{ 49, 6, 5 }, { 49, 6, 25 }, { 49, 9, 5 }, { 49, 9, 25 },
			{ 49, 10, 3 }, { 49, 10, 27 }, { 49, 11, 15 }, { 49, 14, 15 },
			{ 53, 1, 7 }, { 53, 1, 13 }, { 53, 1, 17 }, { 53, 1, 23 },
			{ 53, 4, 7 }, { 53, 4, 13 }, { 53, 4, 17 }, { 53, 4, 23 },
			{ 53, 11, 7 }, { 53, 11, 13 }, { 53, 11, 17 }, { 53, 11, 23 },
			{ 53, 14, 7 }, { 53, 14, 13 }, { 53, 14, 17 }, { 53, 14, 23 } };
	
	private static int[][] segs = new int[60][];
	final static int THRESHOLD = 79000000;
	private static final int LIMIT = 205 * 1000 * 9923;

	private static final int[][] DF2 = { { 7, 1, 2 }, { 7, 1, 8 }, { 7, 1, 22 },
			{ 7, 1, 28 }, { 7, 3, 10 }, { 7, 3, 20 }, { 7, 7, 10 },
			{ 7, 7, 20 }, { 7, 9, 2 }, { 7, 9, 8 }, { 7, 9, 22 }, { 7, 9, 28 },
			{ 19, 1, 4 }, { 19, 1, 14 }, { 19, 1, 16 }, { 19, 1, 26 },
			{ 19, 5, 2 }, { 19, 5, 8 }, { 19, 5, 22 }, { 19, 5, 28 },
			{ 19, 9, 4 }, { 19, 9, 14 }, { 19, 9, 16 }, { 19, 9, 26 },
			{ 31, 3, 2 }, { 31, 3, 8 }, { 31, 3, 22 }, { 31, 3, 28 },
			{ 31, 5, 4 }, { 31, 5, 14 }, { 31, 5, 16 }, { 31, 5, 26 },
			{ 31, 7, 2 }, { 31, 7, 8 }, { 31, 7, 22 }, { 31, 7, 28 },
			{ 43, 1, 10 }, { 43, 1, 20 }, { 43, 3, 4 }, { 43, 3, 14 },
			{ 43, 3, 16 }, { 43, 3, 26 }, { 43, 7, 4 }, { 43, 7, 14 },
			{ 43, 7, 16 }, { 43, 7, 26 }, { 43, 9, 10 }, { 43, 9, 20 } };
	
	private static final int[][] DF3 = { { 11, 0, 7 }, { 11, 0, 13 }, { 11, 0, 17 },
			{ 11, 0, 23 }, { 11, 2, 1 }, { 11, 2, 11 }, { 11, 2, 19 },
			{ 11, 2, 29 }, { 11, 3, 4 }, { 11, 3, 14 }, { 11, 3, 16 },
			{ 11, 3, 26 }, { 11, 5, 2 }, { 11, 5, 8 }, { 11, 5, 22 },
			{ 11, 5, 28 }, { 11, 7, 4 }, { 11, 7, 14 }, { 11, 7, 16 },
			{ 11, 7, 26 }, { 11, 8, 1 }, { 11, 8, 11 }, { 11, 8, 19 },
			{ 11, 8, 29 }, { 23, 1, 10 }, { 23, 1, 20 }, { 23, 2, 7 },
			{ 23, 2, 13 }, { 23, 2, 17 }, { 23, 2, 23 }, { 23, 3, 2 },
			{ 23, 3, 8 }, { 23, 3, 22 }, { 23, 3, 28 }, { 23, 4, 5 },
			{ 23, 4, 25 }, { 23, 6, 5 }, { 23, 6, 25 }, { 23, 7, 2 },
			{ 23, 7, 8 }, { 23, 7, 22 }, { 23, 7, 28 }, { 23, 8, 7 },
			{ 23, 8, 13 }, { 23, 8, 17 }, { 23, 8, 23 }, { 23, 9, 10 },
			{ 23, 9, 20 }, { 47, 1, 4 }, { 47, 1, 14 }, { 47, 1, 16 },
			{ 47, 1, 26 }, { 47, 2, 5 }, { 47, 2, 25 }, { 47, 3, 10 },
			{ 47, 3, 20 }, { 47, 4, 1 }, { 47, 4, 11 }, { 47, 4, 19 },
			{ 47, 4, 29 }, { 47, 6, 1 }, { 47, 6, 11 }, { 47, 6, 19 },
			{ 47, 6, 29 }, { 47, 7, 10 }, { 47, 7, 20 }, { 47, 8, 5 },
			{ 47, 8, 25 }, { 47, 9, 4 }, { 47, 9, 14 }, { 47, 9, 16 },
			{ 47, 9, 26 }, { 59, 0, 1 }, { 59, 0, 11 }, { 59, 0, 19 },
			{ 59, 0, 29 }, { 59, 1, 2 }, { 59, 1, 8 }, { 59, 1, 22 },
			{ 59, 1, 28 }, { 59, 4, 7 }, { 59, 4, 13 }, { 59, 4, 17 },
			{ 59, 4, 23 }, { 59, 5, 4 }, { 59, 5, 14 }, { 59, 5, 16 },
			{ 59, 5, 26 }, { 59, 6, 7 }, { 59, 6, 13 }, { 59, 6, 17 },
			{ 59, 6, 23 }, { 59, 9, 2 }, { 59, 9, 8 }, { 59, 9, 22 },
			{ 59, 9, 28 } };

	private static final int[] dAll = { 1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
			47, 49, 53, 59 };

	private static final int[] Under60 = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41,
			43, 47, 53, 59 };

	private static void enumber2(int f, int g, int d, long L, int B) {
		int x = f, y0 = g;
		long temp = L + B;
		int k0 = (3 * f * f + g * g - d) / 60;

		while (k0 < temp) {
			k0 += x + 5;
			x += 10;
		}

		while (true) {
			x -= 10;
			k0 -= x + 5;
			if (x <= 0)
				break;

			while (k0 < L) {
				k0 += y0 + 15;
				y0 += 30;
			}

			int k = k0, y = y0;
			while (k < temp) {
				segs[d][(int) (k - L) >> 5] ^= 1 << ((k - L) & 31);
				k += y + 15;
				y += 30;
			}
		}
	}

	private static void enumber3(int f, int g, int d, long L, int B) {
		int x = f, y0 = g;
		long temp = L + B;
		int k0 = (3 * f * f - g * g - d) / 60;

		outer: while (true) {
			while (k0 >= temp) {
				if (x <= y0)
					break outer;
				k0 -= y0 + 15;
				y0 += 30;
			}

			int k = k0, y = y0;
			while (k >= L && y < x) {
				segs[d][(int) (k - L) >> 5] ^= 1 << ((k - L) & 31);
				k -= y + 15;
				y += 30;
			}

			k0 += x + 5;
			x += 10;
		}
	}
	// GCD extend
	public static long GCDExtend(long a, long b) {
		long r = 0, s = 1;

		while (b != 0) {
			long c = a / b, d;
			d = a;
			a = b;
			b = d % b;
			d = r;
			r = s;
			s = d - c * s;
		}

		return r;
	}

	// tra ve chi so cua phan tu trong mang 'arr' gan x nhat
	private static int indexingOfNearestFactor(int[] arr, int x) {
		int l = 0, r = arr.length - 1;

		while (l <= r) {
			int c = (l + r) >>> 1;
			int currDiff = arr[c] - x;

			if (currDiff < 0)
				l = c + 1;
			else if (currDiff == 0)
				return c;
			else
				r = c - 1;
		}

		return l;
	}

	public static int[] eratosthenesSieve(int n) {

		int nn = (((n >>> 1) - 1) >>> 5) + 1;
		int[] arr = new int[nn];
		Arrays.fill(arr, 0xffffffff);
		arr[nn - 1] = (1 << (((n >>> 1) - 1) & 31)) - 1;

		int sq = (int) (Math.sqrt(n) - 3) >>> 1;
		for (int p = 0; p <= sq; p++) {
			if ((arr[p >>> 5] & (1 << (p & 31))) != 0) {
				int m = (p << 1) + 3, m2 = m << 1;
				for (int mm = m * m; mm <= n; mm += m2) {
					int ind = (mm - 3) >>> 1;
					arr[ind >>> 5] &= ~(1 << (ind & 31));
				}
			}
		}

		int u = n + 32, pos = 1;
		double lu = Math.log(u);
		int[] ret = new int[(int) (u / lu + u / (lu * lu) * 1.5)];
		ret[0] = 2;

		for (int i = 0; i < nn; i++) {
			for (int j = 0; j <= 31; j++) {
				if ((arr[i] & (1 << j)) != 0)
					ret[pos++] = (((i << 5) | j) << 1) + 3;
			}
		}

		return Arrays.copyOf(ret, pos);
	}

	private static int[] atkinSieve(int n) {
		if (n <= 60) {
			int index = Arrays.binarySearch(Under60, (int) n);
			if (index < 0)
				index = -(index + 2);
			return Arrays.copyOf(Under60, index + 1);
		}

		int B = 60 * (int) Math.sqrt(n);
		int[] primes = eratosthenesSieve((int) Math.sqrt(n));
		int u = n + 32;
		double lu = Math.log(u);
		int[] ret = new int[(int) (u / lu + u / (lu * lu) * 1.5)];
		int r = Under60.length;
		System.arraycopy(Under60, 0, ret, 0, r);

		for (int d : dAll)
			segs[d] = new int[(B >> 5) + 1];

		for (int L = 1; L <= n / 60; L += B) {
			for (int d : dAll)
				Arrays.fill(segs[d], 0);
			int limit = 60 * (L + B);

			for (int[] D : DF1)
				enumber1(D[1], D[2], D[0], L, B); // Case 1
			for (int[] D : DF2)
				enumber2(D[1], D[2], D[0], L, B); // Case 2
			for (int[] D : DF3)
				enumber3(D[1], D[2], D[0], L, B); // Case 3

			for (int p : primes) {
				if ((long) p * p > limit)
					break;
				if (p >= 7) {
					long b = -GCDExtend(p * p, 60);
					int p2 = p * p;
					if (b < 0)
						b += p2;
					for (int d : dAll) {
						int x = (int) (b * (60 * L + d) % p2);
						for (; x < B; x += p2) {
							segs[d][x >> 5] &= ~(1 << (x & 31));
						}
					}
				}
			}

			inner: for (int j = 0; j < (B >> 5) + 1; j++) {
				for (int x = 0; x < 32; x++) {
					long base = 60 * (L + x + (j << 5));
					for (int d : dAll) {
						if (base + d > n)
							break inner;
						if (segs[d][j] << 31 - x < 0) {
							ret[r++] = 60 * (L + x + (j << 5)) + d;
						}
					}
				}
			}
		}
		return Arrays.copyOf(ret, r);
	}

	static {
		int primes[] = SievePrimeNumber.sievePrime(10000), pos = 0;
		smallPrimes = new BigInteger[primes.length];

		for (int i : primes) 
			smallPrimes[pos++] = new BigInteger(Integer.toString(i));
	}


	private static long[] segmentedSieveofEratosthenes(long lo, long hi) {
		int maxPrime = (int) (Math.sqrt(hi)); // Sieve until the square root
		int[] baseSieve = SievePrimeNumber.sievePrime(maxPrime);

		int pos = 0;
		long[] primes = new long[(int) (Math.ceil(1.5 * hi / Math.log(hi)) - Math
				.ceil(1.5 * lo / Math.log(lo)))];

		if (lo < maxPrime) {
			int loPos = indexingOfNearestFactor(baseSieve, (int) lo);
			for (int i = loPos; i < baseSieve.length; i++) {
				primes[pos++] = baseSieve[i];
			}

			lo = maxPrime;
		}

		int length = ((int) (hi - lo) >> 4) + 1;
		if ((lo & 1) == 0)
			lo++;
		byte[] sieve = new byte[length];
		long end = hi;
		hi = lo + (length << 4) - 2;

		int intervalSize = length << 3;
		for (int i = 1; i < baseSieve.length; i++) {
			int p = baseSieve[i];
			int offset = (p - (int) (lo % p)) % p;
			if ((offset & 1) == 1)
				offset += p;

			offset >>= 1;
			for (; offset < intervalSize; offset += p) {
				sieve[offset >> 3] |= (1 << (offset & 7));
			}
		}

		for (long n = lo; n <= end; n += 2) {
			int dn = (int) (n - lo);
			if (((sieve[dn >> 4] >> ((dn >> 1) & 7)) & 1) == 0)
				primes[pos++] = n;
		}

		return Arrays.copyOf(primes, pos);
	}

	private static int[] segmentedSieveofEratosthenes(int lo, int hi) {
		int maxPrime = (int) (Math.sqrt(hi)); 
		int[] baseSieve = SievePrimeNumber.sievePrime(maxPrime);

		int pos = 0;
		int[] primes = new int[(int) (Math.ceil(1.5 * hi / Math.log(hi)) - Math
				.ceil(1.5 * lo / Math.log(lo)))];

		if (lo < maxPrime) {
			int loPos = indexingOfNearestFactor(baseSieve, lo);
			for (int i = loPos; i < baseSieve.length; i++) {
				primes[pos++] = baseSieve[i];
			}

			lo = maxPrime;
		}

		int length = ((hi - lo) >>> 4) + 1;
		if ((lo & 1) == 0)
			lo++;
		byte[] sieve = new byte[length];
		long end = hi;
		hi = lo + (length << 4) - 2;

		int intervalSize = length << 3;
		for (int i = 1; i < baseSieve.length; i++) {
			int p = baseSieve[i];
			int offset = (p - (lo % p)) % p;

			if ((offset & 1) == 1)
				offset += p;

			offset >>= 1;
			for (; offset < intervalSize; offset += p) {
				sieve[offset >> 3] |= (1 << (offset & 7));
			}
		}

		for (int n = lo; n <= end; n += 2) {
			int dn = n - lo;
			if (((sieve[dn >> 4] >> ((dn >> 1) & 7)) & 1) == 0)
				primes[pos++] = n;
		}

		return Arrays.copyOf(primes, pos);
	}


	public static int[] sievePrime(int n) {
		return (n < THRESHOLD) ? eratosthenesSieve(n)
				: atkinSieve(n);
	}

	public static long[] sievePrime(long lo, long hi) {
		if (lo > hi)
			return null;
		if (lo < 2)
			lo = 2;
		return segmentedSieveofEratosthenes(lo, hi);
	}

	public static int[] sievePrime(int lo, int hi) {
		if (lo > hi)
			return null;
		if (lo < 2)
			lo = 2;

		if (hi <= LIMIT) {
			if (hi - lo >= THRESHOLD) {
				int[] primes = atkinSieve(hi);
				return Arrays.copyOfRange(primes,indexingOfNearestFactor(primes, lo), primes.length);
			}
		}
		return segmentedSieveofEratosthenes(lo, hi);
	}

}
