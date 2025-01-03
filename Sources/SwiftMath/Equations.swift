//
//  Equations.swift
//
//  Created by Jacob W Esselstyn on 12/29/22.
//

import SwiftUI

public protocol MathEquation: CustomStringConvertible {
	func callAsFunction(_ x: Double) -> Double
	func makeDerivative() -> MathEquation
	func integrate(plus c: Double) -> MathEquation
}
extension MathEquation {
	
}

public struct Math {
	
	public struct MultiEQ: MathEquation {
		
		public init(segments: [Segment]) {
			self.segments = segments
		}
		public init(equations: [any MathEquation]) {
			var newSegments = [Segment]()
			for eq in equations {
				newSegments.append(Segment(eq: eq))
			}
			self.segments = newSegments
		}
		
		public struct Segment {
			var eq: any MathEquation
			let xStart: Double
			let xEnd: Double
			let xStartIsInclusive: Bool
			let xEndIsInclusive: Bool
			public init(eq: any MathEquation, xStart: Double = Double.leastNonzeroMagnitude, xEnd: Double = Double.greatestFiniteMagnitude, xStartIsInclusive: Bool = true, xEndIsInclusive: Bool = true) {
				self.eq = eq
				self.xStart = xStart
				self.xEnd = xEnd
				self.xStartIsInclusive = xStartIsInclusive
				self.xEndIsInclusive = xEndIsInclusive
			}
		}
		public let segments: [Segment]
		
		public func callAsFunction(_ x: Double) -> Double {
			var result = 0.0
			for s in segments {
				if s.xStart == x && !s.xStartIsInclusive { continue }
				if s.xEnd == x && !s.xEndIsInclusive { continue }
				if s.xStart <= x && s.xEnd >= x { result += s.eq(x) }
			}
			return result
		}
		
		public func makeDerivative() -> any MathEquation {
			var newSegments: [Segment] = []
			for s in segments {
				var segment = s
				segment.eq = s.eq.makeDerivative()
				newSegments.append(segment)
			}
			return MultiEQ(segments: newSegments)
		}
		
		public func integrate(plus c: Double) -> any MathEquation {
			var newSegments: [Segment] = []
			for s in segments {
				var segment = s
				segment.eq = s.eq.integrate(plus: c)
				newSegments.append(segment)
			}
			return MultiEQ(segments: newSegments)
		}
		
		public var description: String {
			var str = "{\n"
			for s in segments {
				str += s.eq.description + ", [\(s.xStart),\(s.xEnd)] \n"
			}
			return str+"\n}"
		}
		
		public func minmaxValues(in range: ClosedRange<Double>?, num_dx: Double = 10) -> (min: Double, max: Double) {
			guard segments.count > 0 else { return (0, 0) }
			let minX = segments.min(by: {$0.xStart < $1.xStart})!.xStart
			let maxX = segments.max(by: {$0.xEnd < $1.xEnd})!.xEnd
			let xStart = max(range?.lowerBound ?? minX, minX)
			let xEnd = min(range?.upperBound ?? maxX, maxX)
			var minValue = callAsFunction(xStart)
			var maxValue = minValue
			for i in 1..<Int(num_dx) {
				let z = Double(i)/num_dx // Percent
				if z > 1 || z < 0 { fatalError("z out of range") }
				let x = xStart + (xEnd - xStart)*z
				minValue = min(minValue, callAsFunction(x))
				maxValue = max(maxValue, callAsFunction(x))
			}
			return (minValue, maxValue)
		}
	}
	
	public struct BasicPolynomialEQ: MathEquation {
		
		public let terms: [Term]
		
		public init(terms: [Term]) {
			self.terms = terms
		}
		public init(terms: Term...) {
			self.terms = terms
		}
//		public init(_ str1: String) {
//			var str2 = str1.trimmingCharacters(in: CharacterSet(arrayLiteral: " "))
//			str2 = str2.replacingOccurrences(of: "-", with: "+-")
//			var termStrs = str2.components(separatedBy: CharacterSet(arrayLiteral: "+"))
//			termStrs.updateEach { s in
//				s.trimmingCharacters(in: CharacterSet(charactersIn: " "))
//			}
//			termStrs.removeAll { $0.isEmpty }
//			var terms = [Term]()
//			for t in termStrs {
//				let t2 = t.replacingOccurrences(of: "x^", with: "x")
//				let components = t2.components(separatedBy: CharacterSet(arrayLiteral: "x"))
//				if components.count == 2 {
//					let a = components.first!
//					let n = components.last!
//					terms.append(Term(Double(a) ?? 1, xToTheDouble: Double(n) ?? 1))
//				} else if components.count == 1 {
//					let a = components.first!
//					terms.append(Term(Double(a) ?? 1, xToThe: 0))
//				}
//			}
//			self.terms = terms
//		}
//		public init(_ polynomial: String) {
//			let pattern = #"([-+]?[0-9]*\.?[0-9]*)x\^?([0-9]*)"#
//			let regex = try! NSRegularExpression(pattern: pattern)
//			let matches = regex.matches(in: polynomial, range: NSRange(polynomial.startIndex..., in: polynomial))
//
//			var terms = [Term]()
//			for match in matches {
//				let coefficient = Double((Range(match.range(at: 1), in: polynomial).flatMap { String(polynomial[$0]) }) ?? "1") ?? 1.0
//				let degree = Double((Range(match.range(at: 2), in: polynomial).flatMap { String(polynomial[$0]) }) ?? "1") ?? 1.0
//				terms.append(Term(coefficient, xToTheDouble: degree))
//			}
//			self.terms = terms
//		}
		
		public func callAsFunction(_ x: Double) -> Double {
			var result = 0.0
			for term in terms {
				result += (term.coefficient*pow(x, term.degree))
			}
			return result
		}
		
		public func makeDerivative() -> MathEquation {
			var newTerms = [Term]()
			for term in terms {
				if term.degree > 0 {
					newTerms.append(Term(term.coefficient*term.degree, xToTheDouble: term.degree-1))
				}
			}
			return BasicPolynomialEQ(terms: newTerms)
		}
		public func integrate(plus c: Double) -> MathEquation {
			var newTerms = [Term]()
			for term in terms {
				if term.degree >= 0 {
					newTerms.append(Term(term.coefficient/(term.degree+1), xToTheDouble: term.degree+1))
				}
			}
			newTerms.append(Term(c, xToTheDouble: 0))
			return BasicPolynomialEQ(terms: newTerms)
		}
		
		public var description: String {
			let termDescriptions = terms.map { term in
				term.description
			}
			return termDescriptions.joined(separator: " + ")
		}
		public struct Term: CustomStringConvertible {
			let degree: Double
			let coefficient: Double
			
			public init(_ coefficient: Double, xToThe degree: UInt) {
				self.degree = Double(degree)
				self.coefficient = coefficient
			}
			init(_ coefficient: Double, xToTheDouble degree: Double) {
				self.degree = Double(degree)
				self.coefficient = coefficient
			}
			
			public var description: String {
				if degree == 1 {
					return "\(coefficient)x"
				} else if degree == 0 {
					return "\(coefficient)"
				} else {
					return "\(coefficient)x^\(degree)"
				}
			}
		}
	}
	
	public struct Impulse: MathEquation {
		
		public let term: Double
		
		public init(term: Double) {
			self.term = term
		}
		
		public func callAsFunction(_ x: Double) -> Double {
			return term
		}
		
		public func makeDerivative() -> MathEquation {
			return Impulse(term: 0)
		}
		public func integrate(plus c: Double) -> MathEquation {
			return BasicPolynomialEQ(terms: BasicPolynomialEQ.Term(term, xToThe: 0))
		}
		
		public var description: String {
			return "\(term) (impulse)"
		}
	}
}
public struct AnyMathEquation: MathEquation {
	private let _callAsFunction: (Double) -> Double
	private let _makeDerivative: () -> MathEquation
	private let _integrate: (Double) -> MathEquation
	private let _description: () -> String

	public init<E: MathEquation>(_ equation: E) {
		_callAsFunction = equation.callAsFunction
		_makeDerivative = equation.makeDerivative
		_integrate = equation.integrate
		_description = { equation.description }
	}

	public func callAsFunction(_ x: Double) -> Double {
		_callAsFunction(x)
	}

	public func makeDerivative() -> MathEquation {
		_makeDerivative()
	}

	public func integrate(plus c: Double) -> MathEquation {
		_integrate(c)
	}

	public var description: String {
		_description()
	}
}
