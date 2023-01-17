double nu = 0.3;

List<Func<double, List<double>, double>> system(int n)
{
    return new List<Func<double, List<double>, double>> {
    (x, y) => -y[1],
    (x, y) => -nu * (n * n / x / x * y[0] + y[1] / x) + y[2],
    (x, y) => -1 / x * ((2 - nu * nu) * n * n / x / x * y[0] + (1 - nu * nu + n * n) / x * y[1] + (1 - nu) * y[2]) + y[3],
    (x, y) => 0
};
}

// See https://aka.ms/new-console-template for more information
Console.WriteLine("Hello, World!");
