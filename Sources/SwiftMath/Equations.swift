//
//  Equations.swift
//
//  Created by Jacob W Esselstyn on 12/29/22.
//

import SwiftUI
import Accelerate

public protocol MathEquation: CustomStringConvertible, Sendable {
	func callAsFunction(_ x: Double) -> Double
	func makeDerivative() -> MathEquation
	func integrate(plus c: Double) -> MathEquation
	mutating func scale(by rhs: Double)
	func scaled(by rhs: Double) -> Self
}
extension MathEquation {
	func definiteIntegral(over domain: ClosedRange<Double>, absoluteTolerance: Double = 1.0e-8, relativeTolerance: Double = 1.0e-2) throws -> Double {
		let q = Quadrature(integrator: .qags(maxIntervals: 12), absoluteTolerance: absoluteTolerance, relativeTolerance: relativeTolerance)
		let result = q.integrate(over: domain, integrand: callAsFunction)
		switch result {
		case .success((let integralResult, let estimatedAbsoluteError)):
			print("quadrature success:", integralResult,
				  estimatedAbsoluteError)
			return integralResult
		case .failure(let error):
			print("quadrature error:", error.errorDescription)
			throw error
		}
	}
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
		
		public struct Segment: Sendable {
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
			public mutating func scale(by rhs: Double) {
				eq.scale(by: rhs)
			}
		}
		public private(set) var segments: [Segment]
		
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
				if s.eq is Math.Impulse { continue }
				var segment = s
				segment.eq = s.eq.makeDerivative()
				newSegments.append(segment)
			}
			return MultiEQ(segments: newSegments)
		}
		
		public func integrate(plus c: Double) -> any MathEquation {
			var newSegments: [Segment] = []
			let sortedSegments = segments.sorted(by: {
				if $0.xStart == $1.xStart {
					return $0.xStartIsInclusive
				} else {
					return $0.xStart<$1.xStart
				}
			})
			// Intergrate segments
			for s in sortedSegments {
				if s.eq is Math.Impulse {
					let segment = Segment(eq: Math.BasicPolynomialEQ(terms: .init(s.eq(0), xToThe: 0)), xStart: s.xStart, xEnd: Double.infinity, xStartIsInclusive: false, xEndIsInclusive: false)
					newSegments.append(segment)
				} else {
					var segment = s
					segment.eq = s.eq.integrate(plus: 0)
					//let cForSegment = MultiEQ(segments: newSegments)(s.xStart) - (segment.eq(s.xStart)-segment.eq(0))
					let cForSegment = -(segment.eq(s.xStart)-segment.eq(0))
					let segmentC = Segment(eq: Math.BasicPolynomialEQ(terms: .init(cForSegment, xToThe: 0)), xStart: s.xStart, xEnd: s.xEnd, xStartIsInclusive: s.xStartIsInclusive, xEndIsInclusive: s.xEndIsInclusive)
					newSegments.append(segment)
					newSegments.append(segmentC)
				}
			}
			// Adjust by c
			if c != 0 {
				let segment = Segment(eq: Math.BasicPolynomialEQ(terms: .init(c, xToThe: 0)), xStart: -Double.infinity, xEnd: Double.infinity, xStartIsInclusive: false, xEndIsInclusive: false)
				newSegments.append(segment)
			}
			return MultiEQ(segments: newSegments)
		}
		public mutating func scale(by rhs: Double) {
			for i in segments.indices {
				segments[i].eq.scale(by: rhs)
			}
		}
		public func scaled(by rhs: Double) -> Self {
			var newSegments: [Segment] = []
			for segment in segments {
				var newSegment = segment
				newSegment.eq.scale(by: rhs)
				newSegments.append(newSegment)
			}
			return Self.init(segments: newSegments)
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
		
		public private(set) var terms: [Term]
		
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
		
		public mutating func scale(by rhs: Double) {
			for i in terms.indices {
				terms[i] = terms[i]*rhs
			}
		}
		public func scaled(by rhs: Double) -> Self {
			var terms: [Term] = []
			for t in self.terms {
				terms.append(t*rhs)
			}
			return Self.init(terms: terms)
		}
		
		public var description: String {
			let termDescriptions = terms.map { term in
				term.description
			}
			return termDescriptions.joined(separator: " + ")
		}
		public struct Term: CustomStringConvertible, Sendable {

			private(set) var degree: Double
			private(set) var coefficient: Double
			
			public init(_ coefficient: Double, xToThe degree: UInt) {
				self.degree = Double(degree)
				self.coefficient = coefficient
			}
			init(_ coefficient: Double, xToTheDouble degree: Double) {
				self.degree = Double(degree)
				self.coefficient = coefficient
			}
			public static func * (_ lhs: Self, _ rhs: Double) -> Self {
				return Term(lhs.coefficient*rhs, xToTheDouble: lhs.degree)
			}
			public static func * (_ lhs: Double, _ rhs: Self) -> Self {
				return Term(lhs*rhs.coefficient, xToTheDouble: rhs.degree)
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
		
		public private(set) var term: Double
		
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
		public mutating func scale(by rhs: Double) {
			term *= rhs
		}
		public func scaled(by rhs: Double) -> Self {
			return Impulse(term: term*rhs)
		}
		
		public var description: String {
			return "\(term) (impulse)"
		}
	}
}
//public struct AnyMathEquation: MathEquation {
//	private let _callAsFunction: (Double) -> Double
//	private let _makeDerivative: () -> MathEquation
//	private let _integrate: (Double) -> MathEquation
//	private let _description: () -> String
//
//	public init<E: MathEquation>(_ equation: E) {
//		_callAsFunction = equation.callAsFunction
//		_makeDerivative = equation.makeDerivative
//		_integrate = equation.integrate
//		_description = { equation.description }
//	}
//
//	public func callAsFunction(_ x: Double) -> Double {
//		_callAsFunction(x)
//	}
//
//	public func makeDerivative() -> MathEquation {
//		_makeDerivative()
//	}
//
//	public func integrate(plus c: Double) -> MathEquation {
//		_integrate(c)
//	}
//
//	public var description: String {
//		_description()
//	}
//}
