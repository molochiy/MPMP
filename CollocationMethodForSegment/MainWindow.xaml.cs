using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows;
using MathNet.Numerics.LinearAlgebra;
using ZedGraph;

namespace CollocationMethodForSegment
{
  /// <summary>
  /// Interaction logic for MainWindow.xaml
  /// </summary>
  public partial class MainWindow : Window
  {
    private InputData _inputData;

    public MainWindow()
    {
      InitializeComponent();
      InitializeInputData();
      InitializeGraphPane();
    }

    private void InitializeInputData()
    {
      GTxbx.Text = "1";
      NTxbx.Text = "10";
    }

    public void InitializeGraphPane()
    {
      FunctionResultGraph.GraphPane.XAxis.Scale.Min = 0;
      FunctionResultGraph.GraphPane.XAxis.Scale.Max = 1;
      FunctionResultGraph.GraphPane.YAxis.Scale.Min = 0;
      FunctionResultGraph.GraphPane.YAxis.Scale.Max = 10;
      FunctionResultGraph.GraphPane.Title.Text = "";
      FunctionResultGraph.GraphPane.XAxis.Title.Text = "t";
      FunctionResultGraph.GraphPane.YAxis.Title.Text = "τ(t)";
      FunctionResultGraph.IsShowPointValues = true;
      FunctionResultGraph.AxisChange();
      FunctionResultGraph.Invalidate();
    }

    private void CalculateBtn_Click(object sender, RoutedEventArgs e)
    {
      InitializeProblem();
      var points = Calculate();

      FillFunctionResultDtGrd(points);
      DrawFunctionResultGraph(points);
    }

    private void InitializeProblem()
    {
      _inputData = new InputData(GTxbx.Text[0] == '-' ? "0" + GTxbx.Text : GTxbx.Text, NTxbx.Text);
    }

    private List<KeyValuePair<double, double>> Calculate()
    {
      var matrix = InitializeMatrix();
      var vector = InitializeVector();

      var result = matrix.Solve(vector);

      var points = new List<KeyValuePair<double, double>>();
      for (int i = 0; i < _inputData.ColocationPoints.Count; i++)
      {
        points.Add(new KeyValuePair<double, double>(_inputData.ColocationPoints[i], result[i]));
      }

      return points;
    }

    private Matrix<double> InitializeMatrix()
    {
      return Matrix<double>.Build.Dense(_inputData.ColocationPoints.Count, _inputData.ColocationPoints.Count,
        (k, j) => MatrixFuntion(k, j + 1));
    }

    private double MatrixFuntion(int k, int j)
    {
      if (_inputData.PartitionPoints[j - 1] > _inputData.ColocationPoints[k])
      {
        return -1 / (2 * Math.PI) * ((_inputData.PartitionPoints[j] - _inputData.ColocationPoints[k]) * Math.Log(_inputData.PartitionPoints[j] - _inputData.ColocationPoints[k]) - (_inputData.PartitionPoints[j - 1] - _inputData.ColocationPoints[k]) * Math.Log(_inputData.PartitionPoints[j - 1] - _inputData.ColocationPoints[k]) - (_inputData.PartitionPoints[j] - _inputData.PartitionPoints[j - 1]));
      }

      if (_inputData.PartitionPoints[j] < _inputData.ColocationPoints[k])
      {
        return -1 / (2 * Math.PI) * (-(-_inputData.PartitionPoints[j] + _inputData.ColocationPoints[k]) * Math.Log(-_inputData.PartitionPoints[j] + _inputData.ColocationPoints[k]) + (-_inputData.PartitionPoints[j - 1] + _inputData.ColocationPoints[k]) * Math.Log(-_inputData.PartitionPoints[j - 1] + _inputData.ColocationPoints[k]) - (_inputData.PartitionPoints[j] - _inputData.PartitionPoints[j - 1]));
      }

      return -1 / (2 * Math.PI) * ((-_inputData.PartitionPoints[j - 1] + _inputData.ColocationPoints[k]) * Math.Log(-_inputData.PartitionPoints[j - 1] + _inputData.ColocationPoints[k]) + (_inputData.PartitionPoints[j] - _inputData.ColocationPoints[k]) * Math.Log(_inputData.PartitionPoints[j] - _inputData.ColocationPoints[k]) - (_inputData.PartitionPoints[j] - _inputData.PartitionPoints[j - 1]));
    }

    private Vector<double> InitializeVector()
    {
      return Vector<double>.Build.Dense(_inputData.ColocationPoints.Count, i => _inputData.G.Calculate(_inputData.ColocationPoints[i]));
    }

    private void FillFunctionResultDtGrd(List<KeyValuePair<double, double>> points)
    {
      var results = points.Select(point => new { T = point.Key, Result = point.Value });
      FunctionResultDtGrd.ItemsSource = results;
    }

    private void DrawFunctionResultGraph(List<KeyValuePair<double, double>> points)
    {
      FunctionResultGraph.GraphPane.CurveList.Clear();
      PointPairList numList = new PointPairList();
      FunctionResultGraph.GraphPane.XAxis.Scale.Min = _inputData.LeftLimit;
      FunctionResultGraph.GraphPane.XAxis.Scale.Max = _inputData.RightLimit;
      double min = points.Min(p => p.Value);
      double max = points.Max(p => p.Value);
      FunctionResultGraph.GraphPane.YAxis.Scale.Min = min - 0.25 * Math.Abs(min);
      FunctionResultGraph.GraphPane.YAxis.Scale.Max = max + 0.25 * Math.Abs(max);
      FunctionResultGraph.AxisChange();
      FunctionResultGraph.Invalidate();
      for (int i = 0; i < points.Count; ++i)
      {
        numList.Add(points[i].Key, points[i].Value);
      }
      numList.Add(points[points.Count - 1].Key, points[points.Count - 1].Value);
      FunctionResultGraph.GraphPane.AddCurve("", numList, System.Drawing.Color.Blue, SymbolType.Star);
      FunctionResultGraph.AxisChange();
      FunctionResultGraph.Invalidate();
    }
  }
}
