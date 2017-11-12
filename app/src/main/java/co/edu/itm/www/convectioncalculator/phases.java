/*******************************************************************************
 * The MIT License (MIT)
 *
 * Convection calculator app
 * Copyright (C) 2017 Luis F. Cardona S.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 ******************************************************************************/

/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 3.0.8
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package co.edu.itm.www.convectioncalculator;

public final class phases {
  public final static phases iphase_liquid = new phases("iphase_liquid");
  public final static phases iphase_supercritical = new phases("iphase_supercritical");
  public final static phases iphase_supercritical_gas = new phases("iphase_supercritical_gas");
  public final static phases iphase_supercritical_liquid = new phases("iphase_supercritical_liquid");
  public final static phases iphase_critical_point = new phases("iphase_critical_point");
  public final static phases iphase_gas = new phases("iphase_gas");
  public final static phases iphase_twophase = new phases("iphase_twophase");
  public final static phases iphase_unknown = new phases("iphase_unknown");
  public final static phases iphase_not_imposed = new phases("iphase_not_imposed");

  public final int swigValue() {
    return swigValue;
  }

  public String toString() {
    return swigName;
  }

  public static phases swigToEnum(int swigValue) {
    if (swigValue < swigValues.length && swigValue >= 0 && swigValues[swigValue].swigValue == swigValue)
      return swigValues[swigValue];
    for (int i = 0; i < swigValues.length; i++)
      if (swigValues[i].swigValue == swigValue)
        return swigValues[i];
    throw new IllegalArgumentException("No enum " + phases.class + " with value " + swigValue);
  }

  private phases(String swigName) {
    this.swigName = swigName;
    this.swigValue = swigNext++;
  }

  private phases(String swigName, int swigValue) {
    this.swigName = swigName;
    this.swigValue = swigValue;
    swigNext = swigValue+1;
  }

  private phases(String swigName, phases swigEnum) {
    this.swigName = swigName;
    this.swigValue = swigEnum.swigValue;
    swigNext = this.swigValue+1;
  }

  private static phases[] swigValues = { iphase_liquid, iphase_supercritical, iphase_supercritical_gas, iphase_supercritical_liquid, iphase_critical_point, iphase_gas, iphase_twophase, iphase_unknown, iphase_not_imposed };
  private static int swigNext = 0;
  private final int swigValue;
  private final String swigName;
}

