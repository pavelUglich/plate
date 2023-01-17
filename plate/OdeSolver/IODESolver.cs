using System;
using System.Collections.Generic;

namespace RungeKutta
{

    /// <summary>
    /// Интерфейс для решения дифференциальных уравнений
    /// </summary>
    /// <typeparam name="T">тип переменной</typeparam>
    /// <typeparam name="U">тип значения функции</typeparam>
    interface IOdEsolver<T, U>
    {
        /// <summary>
        /// правые части дифференциальных уравнений
        /// </summary>
        List<Func<T, List<U>, U>> Equations { get; } 

        /// <summary>
        /// погрешность вычислений
        /// </summary>
        double Epsilon { get; }

        /// <summary>
        /// метод для решения задачи Коши
        /// </summary>
        /// <param name="a">точка, в которой заданы начальные условия</param>
        /// <param name="b">точка, в которой ищется решение</param>
        /// <param name="initialConditions">вектор начальных условий</param>
        /// <returns>решение задачи Коши в точке x=b</returns>
        List<T> Solve(T a, T b, List<T> initialConditions);
    }
}