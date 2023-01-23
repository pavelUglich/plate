namespace RungeKutta
{
    public class OdeSolver : IOdEsolver<double, double>
    {
        public List<Func<double, List<double>, double>> Equations { get; }
        public double Epsilon { get; }
        private readonly double[][] _butcherTableau;

        /// <summary>
        /// конструктор класса
        /// </summary>
        /// <param name="funcs"> набор правых частей системы уравнений </param>
        /// <param name="butcherTableau"> таблица Бутчера </param>
        /// <param name="epsilon"> погрешность вычислений</param>
        public OdeSolver(List<Func<double, List<double>, double>> funcs,
            double[][] butcherTableau, double epsilon = 0.1e-5)
        {
            Equations = funcs;
            _butcherTableau = butcherTableau;
            Epsilon = epsilon;
        }

        /// <summary>
        /// Метод, строящий решение системы уравнений и возвращающий его в виде
        /// списка вещественных чисел
        /// </summary>
        /// <param name="a"> точка, в которой задаются начальный условия
        /// </param>
        /// <param name="b"> точка, вкоторой ищутся решения </param>
        /// <param name="initialConditions"> начальные условия </param>
        /// <returns> решение задачи Коши </returns>
        public List<double> Solve(double a, double b,
            List<double> initialConditions)
        {
            double h = b - a;
            while (a < b && Math.Abs(b - a) > double.Epsilon)
            {
                var r = CalculateTheResidual(a, initialConditions, h,
                    out List<double> u);
                while (Math.Abs(r) > Epsilon)
                {
                    h /= 2;
                    r = CalculateTheResidual(a, initialConditions, h, out u);
                }
                a += h;
                initialConditions = u;
            }
            return initialConditions;
        }

        /// <summary>
        /// Метод ля вычисления невязки по формуле ()
        /// </summary>
        /// <param name="a"> точка, в которой задаются начальный условия
        /// </param>
        /// <param name="initialConditions"> начальные условия </param>
        /// <param name="h"> шаг метода </param>
        /// <param name="u"> решение в точке a + h </param>
        /// <returns> невязка </returns>
        private double CalculateTheResidual(double a,
            List<double> initialConditions, double h, out List<double> u)
        {
            EvaluateUh(initialConditions, a, h, out u, out var u_);
            List<double> difference = new List<double>(u);
            for (int i = 0; i < u.Count; i++)
            {
                difference[i] -= u_[i];
            }
            return Math.Sqrt(difference.Select(x => x * x).Sum());
        }

        /// <summary>
        /// Находит приближенные решения u и uh и возвращает их в виде выходных
        /// параметров
        /// </summary>
        /// <param name="initialConditions"> начальные условия </param>
        /// <param name="a"> точка, в которой задаются начальный условия
        /// </param>
        /// <param name="h"> шаг метода </param>
        /// <param name="u"> приближенное решение в точке a + h </param>
        /// <param name="uh"> приближенное решение в точке a + h </param>
        private void EvaluateUh(List<double> initialConditions, double a,
            double h, out List<double> u, out List<double> uh)
        {
            u = new List<double>(initialConditions);
            uh = new List<double>(initialConditions);
            var k = EvaluateK(initialConditions, a, h);
            int size = _butcherTableau.Length - 1;
            for (int i = 0; i < initialConditions.Count; i++)
            {
                for (int j = 1; j < _butcherTableau[size].Length; j++)
                {
                    u[i] += _butcherTableau[size - 1][j] * k[j - 1][i];
                    uh[i] += _butcherTableau[size][j] * k[j - 1][i];
                }
            }
        }

        /// <summary>
        /// вычисление коэффициентов K_ij
        /// </summary>
        /// <param name="initialConditions"> начальные условия </param>
        /// <param name="a"> левая граница отрезка </param>
        /// <param name="h"> шаг </param>
        /// <returns> матрица K </returns>
        private List<List<double>> EvaluateK(List<double> initialConditions,
            double a, double h)
        {
            var result = new List<List<double>>
            {
                Equations.Select(x => h * x(a, initialConditions)).ToList()
            };
            for (int i = 0; i < _butcherTableau.Length - 2; i++)
            {
                double xx = a + _butcherTableau[i][0] * h;
                List<double> uh = new List<double>(initialConditions);
                for (int j = 0; j < uh.Count; j++)
                {
                    for (int k = 0; k < result.Count; k++)
                    {
                        uh[j] += _butcherTableau[i][k + 1] * result[k][j];
                    }
                }
                result.Add(Equations.Select(x => h * x(xx, uh)).ToList());
            }
            return result;
        }
    }
}
