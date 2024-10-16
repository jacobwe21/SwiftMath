//
//  Math.swift
//  
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
	}
	
	public struct BasicPolynomialEQ: MathEquation {
		
		public let terms: [Term]
		
		public init(terms: [Term]) {
			self.terms = terms
		}
		public init(terms: Term...) {
			self.terms = terms
		}
		public init(_ str1: String) {
			var str2 = str1.trimmingCharacters(in: CharacterSet(arrayLiteral: " "))
			str2 = str2.replacingOccurrences(of: "-", with: "+-")
			var termStrs = str2.components(separatedBy: CharacterSet(arrayLiteral: "+"))
			termStrs.updateEach { s in
				s.trimmingCharacters(in: CharacterSet(charactersIn: " "))
			}
			termStrs.removeAll { $0.isEmpty }
			var terms = [Term]()
			for t in termStrs {
				let t2 = t.replacingOccurrences(of: "x^", with: "x")
				let components = t2.components(separatedBy: CharacterSet(arrayLiteral: "x"))
				if components.count == 2 {
					let a = components.first!
					let n = components.last!
					terms.append(Term(Double(a) ?? 1, xToTheDouble: Double(n) ?? 1))
				} else if components.count == 1 {
					let a = components.first!
					terms.append(Term(Double(a) ?? 1, xToThe: 0))
				}
			}
			self.terms = terms
		}
		
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
}
