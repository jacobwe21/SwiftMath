// swift-tools-version: 5.10
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "SwiftMath",
	defaultLocalization: "en",
	platforms: [.iOS(.v16),.macOS(.v13),.macCatalyst(.v16)],
    products: [
        // Products define the executables and libraries a package produces, making them visible to other packages.
        .library(
            name: "SwiftMath",
            targets: ["SwiftMath"]),
    ],
	dependencies: [
		// Dependencies declare other packages that this package depends on.
		// .package(url: /* package url */, from: "1.0.0"),
		.package(url: "https://github.com/jacobwe21/MySwift", branch: "main")
	],
    targets: [
        // Targets are the basic building blocks of a package, defining a module or a test suite.
        // Targets can depend on other targets in this package and products from dependencies.
        .target(
            name: "SwiftMath", dependencies: ["MySwift"]),
        .testTarget(
            name: "SwiftMathTests",
            dependencies: ["SwiftMath","MySwift"]),
    ]
)
