using System;
using System.Collections.Generic;
using PoohMathParser;

namespace GalerkinMethodForSegment
{
  class InputData
  {
    public InputData(string g, string n)
    {
      G = new MathExpression(g);
      N = Convert.ToInt32(n);

      H = (RightLimit - LeftLimit) / N;

      InitPartitionPoints();
      InitColocationPoints();
    }

    public int N { get; }

    public MathExpression G { get; }

    public double LeftLimit => -1.0;
    public double RightLimit => 1.0;

    public double H { get; }

    public List<double> PartitionPoints { get; set; }

    public List<double> ColocationPoints { get; set; }

    private void InitPartitionPoints()
    {
      PartitionPoints = new List<double>();

      for (var i = 0; i < N; i++)
      {
        PartitionPoints.Add(LeftLimit + i * H);
      }

      PartitionPoints.Add(RightLimit);
    }

    private void InitColocationPoints()
    {
      ColocationPoints = new List<double>();

      for (var i = 0; i < PartitionPoints.Count - 1; i++)
      {
        ColocationPoints.Add((PartitionPoints[i] + PartitionPoints[i + 1]) / 2.0);
      }
    }
  }
}
