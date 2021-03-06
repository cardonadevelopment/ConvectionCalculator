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

public final class composition_types {
  public final static composition_types IFRAC_MASS = new composition_types("IFRAC_MASS");
  public final static composition_types IFRAC_MOLE = new composition_types("IFRAC_MOLE");
  public final static composition_types IFRAC_VOLUME = new composition_types("IFRAC_VOLUME");
  public final static composition_types IFRAC_UNDEFINED = new composition_types("IFRAC_UNDEFINED");
  public final static composition_types IFRAC_PURE = new composition_types("IFRAC_PURE");

  public final int swigValue() {
    return swigValue;
  }

  public String toString() {
    return swigName;
  }

  public static composition_types swigToEnum(int swigValue) {
    if (swigValue < swigValues.length && swigValue >= 0 && swigValues[swigValue].swigValue == swigValue)
      return swigValues[swigValue];
    for (int i = 0; i < swigValues.length; i++)
      if (swigValues[i].swigValue == swigValue)
        return swigValues[i];
    throw new IllegalArgumentException("No enum " + composition_types.class + " with value " + swigValue);
  }

  private composition_types(String swigName) {
    this.swigName = swigName;
    this.swigValue = swigNext++;
  }

  private composition_types(String swigName, int swigValue) {
    this.swigName = swigName;
    this.swigValue = swigValue;
    swigNext = swigValue+1;
  }

  private composition_types(String swigName, composition_types swigEnum) {
    this.swigName = swigName;
    this.swigValue = swigEnum.swigValue;
    swigNext = this.swigValue+1;
  }

  private static composition_types[] swigValues = { IFRAC_MASS, IFRAC_MOLE, IFRAC_VOLUME, IFRAC_UNDEFINED, IFRAC_PURE };
  private static int swigNext = 0;
  private final int swigValue;
  private final String swigName;
}

