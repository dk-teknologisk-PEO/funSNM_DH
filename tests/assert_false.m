function assert_false(condition, name)
%ASSERT_FALSE Asserts that condition is false.
    if condition
        error('ASSERT_FALSE failed: %s', name);
    end
end