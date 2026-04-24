function assert_near(actual, expected, tol, name)
%ASSERT_NEAR Asserts that actual is within tol of expected.
    if abs(actual - expected) > tol
        error('ASSERT_NEAR failed for %s: expected %.6f, got %.6f (tol=%.6f)', ...
            name, expected, actual, tol);
    end
end