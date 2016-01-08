package ecm;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.math.BigInteger;
import java.security.SecureRandom;

public class EllipticCurvesMethod {

    public class Montgomery {

        public BigInteger x;
        public BigInteger z;

        public Montgomery(BigInteger x, BigInteger z) {
            this.x = x;
            this.z = z;
        }

        public Montgomery() {
            this(BigInteger.ZERO, BigInteger.ZERO);
        }

        public Montgomery(Montgomery m) {
            this.x = m.x;
            this.z = m.z;
        }

        public String toString() {
            return "(" + this.x + ", " + this.z + ")";
        }

    }

    public Montgomery doublePoint(Montgomery p) {
        BigInteger u = p.x.add(p.z);
        BigInteger v = p.x.subtract(p.z);
        BigInteger u2 = u.multiply(u);
        BigInteger v2 = v.multiply(v);
        BigInteger t = u2.subtract(v2);
        BigInteger x = u2.multiply(v2).mod(n);
        BigInteger z = t.multiply(v2.add(A24.multiply(t))).mod(n);
        return new Montgomery(x, z);
    }

    private final double[] v = {1.61803398875, 1.72360679775, 1.618347119656,
        1.617914406529, 1.612429949509, 1.632839806089, 1.620181980807,
        1.580178728295, 1.617214616534, 1.38196601125};

    public Montgomery addAPoint(Montgomery p,
            Montgomery q, Montgomery r) {
        /* u = (x2 - z2)(z1 + z1) */
        BigInteger u = p.x.subtract(p.z).multiply(q.x.add(q.z));
        /* v = (x2 + z2)(x1 - z1) */
        BigInteger v = p.x.add(p.z).multiply(q.x.subtract(q.z));
        BigInteger uPlusV = u.add(v);
        BigInteger uMinusV = u.subtract(v);
        BigInteger x = r.z.multiply(uPlusV.multiply(uPlusV)).mod(n);
        BigInteger z = r.x.multiply(uMinusV.multiply(uMinusV)).mod(n);
        return new Montgomery(x, z);
    }

    public Montgomery multiplyScalar(BigInteger k,
            Montgomery P) {

        String kString = k.toString(2);
        int lenKString = kString.length();
        Montgomery Q = P;
        Montgomery R = doublePoint(P);

        for (int i = 1; i < lenKString; i++) {
            if (kString.charAt(i) == '1') {
                Q = addAPoint(R, Q, P);
                R = doublePoint(R);
            } else {
                R = addAPoint(Q, R, P);
                Q = doublePoint(Q);
            }
        }

        return Q;
    }

    private final static BigInteger ZERO = new BigInteger("0");
    private final static BigInteger ONE = new BigInteger("1");
    private final static BigInteger TWO = new BigInteger("2");
    private final static BigInteger THREE = new BigInteger("3");
    private final static BigInteger FIVE = new BigInteger("5");

    public BigInteger factorize(BigInteger n) {
        if (n.isProbablePrime(40)) {
            return ONE;
        }

        this.n = n;
        BigInteger g = ONE;

        double _q = Math.pow(10, (n.toString().length() - 2) >>> 1);
        double logQ = Math.log(_q);

        int B1 = (int) Math.ceil(Math.exp(Math.sqrt(0.5 * logQ * Math.log(logQ))) / 10) * 10;
        int B2 = B1 * 100;
        int D = (int) Math.sqrt(B2);

        double logB1 = Math.log10(B1);
        int[] primes = SievePrimeNumber.sievePrime(B1);

        BigInteger k = ONE;

        for (int prime : primes) {
            k = k.multiply(BigInteger.valueOf((long) Math.pow(prime, Math.floor(logB1 / Math.log10(prime)))));
        }

        primes = SievePrimeNumber.sievePrime(B1 + 1, B2);
        int numPrimes = primes.length;
        Montgomery[] S = new Montgomery[D + 1];
        BigInteger[] beta = new BigInteger[D + 1];

        int curves = 0;

        while (g.equals(ONE) || g.equals(n)) {
            BigInteger sigma = new BigInteger(32, new SecureRandom());

            BigInteger u = sigma.pow(2).subtract(FIVE).mod(n);
            BigInteger v = sigma.shiftLeft(2).mod(n);
            BigInteger A = v.subtract(u).pow(3)
                    .multiply(THREE.multiply(u).add(v))
                    .divide(u.pow(3).shiftLeft(2).multiply(v)).subtract(TWO)
                    .mod(n);
            A24 = A.add(TWO).shiftRight(2);

            Montgomery P = new Montgomery(u.pow(3)
                    .divide(v.pow(3)).mod(n), ONE);

            System.out.println("Curve " + (++curves) + ": " + sigma);

            Montgomery Q = multiplyScalar(k, P);
            g = n.gcd(Q.z);

            if ((g.equals(ONE)) || (g.equals(n))) {
                S[1] = doublePoint(Q);
                S[2] = doublePoint(S[1]);

                for (int d = 1; d <= D; d++) {

                    if (d > 2) {
                        S[d] = addAPoint(S[d - 1], S[1], S[d - 2]);
                    }

                    beta[d] = S[d].x.multiply(S[d].z).mod(n);
                }

                g = ONE;
                int B = B1 - 1;
                Montgomery R = multiplyScalar(BigInteger.valueOf(B),
                        Q);
                Montgomery T = multiplyScalar(
                        BigInteger.valueOf(B - (D << 1)), Q);

                int q = 0;
                int step = D << 1;
                for (int r = B; r < B2; r += step) {
                    BigInteger alpha = R.x.multiply(R.z).mod(n);
                    int limit = r + step;
                    while ((q < numPrimes) && (primes[q] <= limit)) {
                        int delta = (primes[q] - r) >>> 1;
                        g = g.multiply(
                                R.x.subtract(S[delta].x)
                                .multiply(R.z.add(S[delta].z))
                                .subtract(alpha).add(beta[delta])).mod(
                                        n);
                        q++;
                    }

                    Montgomery tmpR = R;
                    R = addAPoint(R, S[D], T);
                    T = tmpR;
                }

                g = n.gcd(g);
            }
        }

        return g;
    }

    private BigInteger n;
    private BigInteger A24;

    private final static int costADD = 6;

    private final static int costDouble = 5;

    public static void main(String args[]) throws IOException {
        EllipticCurvesMethod obj = new EllipticCurvesMethod();

        BigInteger num = new BigInteger("459344432904829185748347139878533089");
        System.out.println("So nguyen N = " + num);
        long time = System.nanoTime();

        boolean P = num.isProbablePrime(1);
        int it = 1;
        while (P == false) {
            BigInteger div = obj.factorize(num);
            System.out.println("Thua so nguyen to " + it + " :" + div);
            num = num.divide(div);
            P = num.isProbablePrime(1);
            it++;
        }
        System.out.println("Thua so nguyen to " + it + ": " + num);

    }
}
