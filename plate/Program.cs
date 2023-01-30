using bem;
using bem.BoundaryValueProblem;
using RungeKutta;

namespace plate
{
    static class Plate
    {
        static double nu = 0.3;
        //static double kappa = 20.0;
        static double innerRadius = 0.00001;

        static List<Func<double, List<double>, double>> system(int n, double kappa)
        {
            return new List<Func<double, List<double>, double>> {
    (x, y) => -y[1],
    (x, y) => -nu*y[1]/x+y[2],
    (x, y) => (1-nu*nu)*y[1]/x/x-(1-nu)*y[2]/x+y[3],
    (x, y) => -y[3]/x - kappa*kappa*y[0]
    };
        }

        /// <summary>
        /// амплитудно-частотная характеристика
        /// </summary>
        /// <param name="maxKappa">наибольшая частота (построение происходит с нуля)</param>
        /// <param name="numPoints">количество узловых точек</param>
        /// <param name="func">функция (АЧХ можно задать аналитическим выражением)</param>
        /// <returns>АЧХ в виде словаря</returns>
        static Dictionary<double, double> FrequencyResponse(double maxKappa,
            int numPoints, Func<double, double> func)
        {
            Dictionary<double, double> result = new Dictionary<double, double>();
            var h = maxKappa / numPoints;
            for (int i = 0; i < numPoints; i++)
            {
                var kappa = i * h;
                result.Add(kappa, func(kappa));
            }
            return result;
        }

        /// <summary>
        /// перегрузка FrequenceResponse с минимальной частотой
        /// </summary>
        /// <param name="minKappa"></param>
        /// <param name="maxKappa"></param>
        /// <param name="numPoints"></param>
        /// <param name="func"></param>
        /// <returns></returns>    
        static Dictionary<double, double> FrequencyResponse(double minKappa,
            double maxKappa, int numPoints, Func<double, double> func)
        {
            Dictionary<double, double> result = new Dictionary<double, double>();
            var h = (maxKappa - minKappa) / numPoints;
            for (int i = 0; i < numPoints; i++)
            {
                var kappa = minKappa + i * h;
                result.Add(kappa, func(kappa));
            }
            return result;
        }

        /// <summary>
        /// частотное уравнение
        /// </summary>
        /// <param name="kappa">максимальная частота</param>
        /// <returns>ачх в виде списка вещественных чисел</returns>
        static double FrecquencyEquation(double kappa)
        {
            var equations = system(0, kappa);
            OdeSolver odeSolver = new OdeSolver(equations, ButcherTableau.RungeKuttaFeldberg);
            var solution1 = odeSolver.Solve(innerRadius, 1.0, new List<double> { 0, 0, 1, 0 });
            var solution2 = odeSolver.Solve(innerRadius, 1.0, new List<double> { 0, 0, 0, 1 });
            Console.WriteLine(kappa);
            return solution1[2] * solution2[3] - solution1[3] * solution2[2];
        }


        /// <summary>
        /// амплитудно-частотная характе
        /// </summary>
        /// <param name="kappa">максимальная частота</param>
        /// <returns>ачх в виде списка вещественных чисел</returns>
        static double FrecquencyResponse(double kappa)
        {
            var equations = system(0, kappa);
            Dictionary<int, double> left =
                new Dictionary<int, double> { { 0, 1 }, { 1, 0 } };
            Dictionary<int, double> right
                = new Dictionary<int, double> { { 2, 0 }, { 3, 0 } };
            BoundaryValueProblem boundaryValueProblem = new BoundaryValueProblem(
                equations, left, right,
                new Dictionary<int, Func<double, double>>(), innerRadius);
            Console.WriteLine(kappa);
            return boundaryValueProblem.Solve()[0];
        }

        /// <summary>
        /// метод секущих
        /// </summary>
        /// <param name="a">левый конец отрезка</param>
        /// <param name="b">правый конец орезка</param>
        /// <param name="equation">уранение</param>
        /// <param name="epsilon">погрешность</param>
        /// <returns>набор корней</returns>
        static double SecantMethod(double a, double b,
            Func<double, double> equation, double epsilon = 0.1e-6)
        {
            while (Math.Abs(a - b) > epsilon)
            {
                var fb = equation(b);
                var c = b - fb * (b - a) / (fb - equation(a));
                (a, b) = (b, c);
            }
            return b;
        }

        /// <summary>
        /// отыскание корней    
        /// </summary>
        /// <param name="maxKappa">максимальная частота</param>
        /// <param name="h">длина отрезка</param>
        /// <param name="equation">уравнение</param>
        /// <returns>набор корней</returns>
        static List<double> Roots(double maxKappa, double h,
            Func<double, double> equation)
        {
            List<double> result = new List<double>();
            double a = 0;
            double b = h;
            double fa = equation(a);
            while (a < maxKappa)
            {
                double fb = equation(b);
                if (fa * fb < 0)
                {
                    result.Add(SecantMethod(a, b, equation));
                }
                a = b;
                fa = fb;
                b += h;
            }
            return result;
        }

        static Dictionary<double, double> EigenMode(double kappa, IEnumerable<double> points)
        {
            var equations = system(0, kappa);
            OdeSolver odeSolver = new OdeSolver(equations, ButcherTableau.RungeKuttaFeldberg);
            var solution1 = odeSolver.Solve(innerRadius, 1.0, new List<double> { 0, 0, 1, 0 });
            var solution2 = odeSolver.Solve(innerRadius, 1.0, new List<double> { 0, 0, 0, 1 });
            List<double> initials = new List<double> { 0, 0, solution2[3], -solution1[3] };
            var solution = odeSolver.Solve(points, initials);
            var maximum = solution.Select(x => x.Value.First()).OrderBy(x => Math.Abs(x)).Last();
            return solution.ToDictionary(x => x.Key, x => x.Value.First() / maximum);
        }


        static void PlotEigenModes(double maxKappa, string fileName, int points = 20)
        {
            var eigenFrecquencies = Roots(maxKappa, 0.5, x => FrecquencyEquation(x));
            Tikz.FileHeader(fileName, out StreamWriter streamWriter);
            var hh = (1 - innerRadius) / points;
            List<double> p = new List<double>();
            for (int i = 0; i <= points; i++)
            {
                p.Add(innerRadius + i * hh);
            }
            streamWriter.WriteLine(@"\legend{");
            foreach (var frecquency in eigenFrecquencies)
            {
                streamWriter.WriteLine($"{frecquency},");
            }
            streamWriter.WriteLine(@"};");
            foreach (var frecquency in eigenFrecquencies)
            {
                var mode = EigenMode(frecquency, p);
                Tikz.AddTheCurve(mode, streamWriter);
            }
            Tikz.FileFooter(streamWriter);
            streamWriter.Close();
        }

        static void Main()
        {
            PlotEigenModes(100, "modes.txt", 50);
        }
    }
}