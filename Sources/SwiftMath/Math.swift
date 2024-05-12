//
//  Math.swift
//  
//
//  Created by Jacob W Esselstyn on 12/29/22.
//

import SwiftUI
import MySwift

public protocol MathEquation: CustomStringConvertible {
	func callAsFunction(_ x: Double) -> Double
	func makeDerivative() -> MathEquation
	func integrate(plus c: Double) -> MathEquation
}
extension MathEquation {
	
}

public struct Math {
	
	public struct PolynomialEQ: MathEquation {
		
		let terms: [Term]
		
		init(terms: [Term]) {
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
					terms.append(Term(Double(a) ?? 1, xToThe: Double(n) ?? 1))
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
				result += (term.coefficient*x**term.degree)
			}
			return result
		}
		
		public func makeDerivative() -> MathEquation {
			var newTerms = [Term]()
			for term in terms {
				if term.degree > 0 {
					newTerms.append(Term(term.coefficient*term.degree, xToThe: term.degree-1))
				}
			}
			return PolynomialEQ(terms: newTerms)
		}
		public func integrate(plus c: Double) -> MathEquation {
			var newTerms = [Term]()
			for term in terms {
				if term.degree > 0 {
					newTerms.append(Term(term.coefficient/(term.degree+1), xToThe: term.degree+1))
				}
			}
			newTerms.append(Term(c, xToThe: 0))
			return PolynomialEQ(terms: newTerms)
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
			
			public init(_ coefficient: Double, xToThe degree: Double) {
				self.degree = degree
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
