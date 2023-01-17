﻿namespace RungeKutta
{
    public static class ButcherTableau
    {
        /// <summary>
        /// матрица Батчера для метода Хойна
        /// </summary>
        public static readonly double[][] Heun = {
            new[]{ 1.0, 1.0 },
            new[]{ 0.0, 0.5, 0.5 },
            new[]{ 0.0, 1.0, 0.0 }
        };

        /// <summary>
        /// Матрица Батчера для метода Богацки-Шампайна
        /// </summary>
        public static readonly double[][] BogackiShampine = {
            new[]{ 0.5,  0.5 },
            new[]{ 0.75, 0.0,      0.75 },
            new[]{ 1.0,  2.0 / 9,  1.0 / 3, 4.0 / 9 },
            new[]{ 0.0,  2.0 / 9,  1.0 / 3, 4.0 / 9, 0.0 },
            new[]{ 0.0,  7.0 / 24, 0.25,    1.0 / 3, 0.125 }
        };

        /// <summary>
        /// Матрица Батчера для метода Кеша-Карпа
        /// </summary>
        public static readonly double[][] CashCarp = {
            new[]{ 0.2, 0.2 },
            new[]{ 0.3, 3.0 / 40,   9.0 / 40 },
            new[]{ 0.6, 0.3,        -0.9,      1.2 },
            new[]{ 1.0, -11.0 / 54, 2.5, -70.0 / 27, 35.0 / 27 },
            new[]{ 0.875, 1631.0 / 55296, 175.0 / 512, 575.0 / 13824,
                44275.0 / 110592, 253.0 / 4096},
            new[]{ 0.0, 37.0 / 378, 0.0, 250.0 / 621, 125.0 / 594, 0.0,
                512.0 / 1771},
            new[]{ 0.0, 2825.0 / 27648, 0.0, 18575.0 / 48384, 13525.0 / 55296,
                277.0 / 14336, 0.25}
        };

        /// <summary>
        /// Матрица Батчера для метода Рунге-Кутты-Фельдберга
        /// </summary>
        public static readonly double[][] RungeKuttaFeldberg = {
            new[]{ 0.25,0.25 },
            new[]{ 0.375, 0.09375, 0.28125 },
            new[]{ 12 / 13.0, 1932 / 2197.0, -7200 / 2197.0, 7296 / 2197.0 },
            new[]{ 1, 439 / 216.0, -8, 3680 / 513.0, -845 / 4104.0 },
            new[]{ 0.5, -8 / 27.0, 2, -3544 / 2565.0, 1859 / 4104.0, -11 / 40.0 },
            new[]{ 0, 16 / 135.0, 0, 6656 / 12825.0, 28561 / 56430.0,
                -9 / 50.0, 2 / 55.0 },
            new[]{ 0, 25 / 216.0, 0, 1408 / 2565.0, 2197 / 4104.0, -0.2, 0 }
        };

        /// <summary>
        /// Матрица Батчера для метода Дорманда-Принса
        /// </summary>
        public static readonly double[][] DormandPrince =
        {
            new[] {0.2, 0.2},
            new[] {0.3, 0.075, 0.225},
            new[] {0.8, 44.0 / 45, -56.0 / 15, 32.0 / 9},
            new[] {8.0 / 9, 19372.0 / 6561, -25360.0 / 2187, 64448.0 / 6561,
                -212.0 / 729},
            new[] {1.0, 9017.0 / 3168, -355.0 / 33, 46732.0 / 5247, 49.0 / 176,
                -5103.0 / 18656},
            new[] {1.0, 35.0 / 384, 0.0, 500.0 / 1113, 125.0 / 192,
                -2187.0 / 6784, 11.0 / 84},
            new[] {0.0, 35.0 / 384, 0.0, 500.0 / 1113, 125.0   / 192,
                -2187.0 / 6784, 11.0 / 84, 0.0},
            new[] {0.0, 5179.0 / 57600, 0.0, 7571.0 / 16695, 393.0 / 640,
                -92097.0 / 339200, 187.0 / 2100, 0.025}
        };

        /// <summary>
        /// Матрица Батчера для метода Дорманда-Принса порядка 7(8)
        /// </summary>
        public static readonly double[][] Rkf78 =
        {
            new[] {1.0 / 18, 1.0 / 18},
            new[] {1.0 / 12, 1.0 / 48, 1.0 / 16},
            new[] {0.125, 0.03125, 0.0, 0.09375},
            new[] {0.3125, 0.3125, 0.0, -1.171875, 1.171875},
            new[] {0.375, 0.0375, 0.0, 0.0, 0.1875, 0.15},
            new[] {0.1475, 29443841.0 / 614563906, 0.0, 0.0,
                77736538.0 / 692538347, -28693883.0 / 1125000000,
                23124283.0 / 1800000000},
            new[] {93.0 / 200, 16016141.0 / 946692911, 0.0, 0.0,
                61564180.0 / 158732637, 22789713.0 / 633445777,
                545815736.0 / 2771057229, -180193667.0 / 1043307555},
            new[] {5490023248.0 / 9719169821, 39632708.0 / 573591083, 0.0, 0.0,
                -433636366.0 / 683701615, -421739975.0 / 2616292301,
                100302831.0 / 723423059, 790204164.0 / 839813087,
                800635310.0 / 3783071287},
            new[] {13.0 / 20, 246121993.0 / 1340847787, 0.0, 0.0,
                -37695042795.0 / 15268766246, -309121744.0 / 1061227803,
                -12992083.0 / 490766935, 6005943493.0 / 2108947869,
                393006217.0 / 1396673457,    123872331.0 / 1001029789},
            new[] {1201146811.0 / 1299019798, -1028468189.0 / 846180014, 0.0,
                0.0, 8478235783.0 / 508512852, 1311729495.0 / 1432422823,
                -10304129995.0 / 1701304382, -48777925059.0 / 3047939560,
                15336726248.0 / 1032824649, -45442868181.0 / 3398467696,
                3065993473.0 / 597172653},
            new[] {1.0, 185892177.0 / 718116043, 0.0, 0.0,
                -3185094517.0 / 667107341, -477755414.0 / 1098053517,
                -703635378.0 / 230739211, 5731566787.0 / 1027545527,
                5232866602.0 / 850066563, -4093664535.0 / 808688257,
                3962137247.0 / 1805957418, 65686358.0 / 487910083},
            new[] {1.0, 403863854.0 / 491063109, 0.0, 0.0,
                -5068492393.0 / 434740067, -411421997.0 / 543043805,
                652783627.0 / 914296604, 11173962825.0 / 925320556,
                -13158990841.0 / 6184727034, 3936647629.0 / 1978049680,
                -160528059.0 / 685178525, 248638103.0 / 1413531060, 0},
            new[] {0.0, 14005451.0 / 335480064, 0.0, 0.0, 0.0, 0.0,
                -59238493.0 / 1068277825, 181606767.0 / 758867731,
                561292985.0 / 797845732, -1041891430.0 / 1371343529,
                760417239.0 / 1151165299, 118820643.0 / 751138087,
                -528747749.0 / 2220607170, 0.25},
            new[] {0.0, 13451932.0 / 455176623, 0.0, 0.0, 0.0, 0.0,
                -808719846.0 / 976000145, 1757004468.0 / 5645159321,
                656045339.0 / 265891186, -3867574721.0 / 1518517206,
                465885868.0 / 322736535, 53011238.0 / 667516719, 2.0 / 45, 0.0}
        };
    }
}