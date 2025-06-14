import XCTest
import MySwiftMath

final class SwiftMathTests: XCTestCase {
    func testExample() throws {
        // XCTest Documentation
        // https://developer.apple.com/documentation/xctest

        // Defining Test Cases and Test Methods
        // https://developer.apple.com/documentation/xctest/defining_test_cases_and_test_methods
    }
}


final class IntegrationTests: XCTestCase {
	func testIntegration() throws {
		let eq = Math.PolynomialEQ(terms: [.init(2,exp:2),.init(1,exp:1),.init(-4,exp:0)])
		let result = try eq.definiteIntegral(over: -6...8)
		XCTAssertEqual(result, 443.3333333333333333333333333333, accuracy: 1e-12)
	}
}
