function assert_true(condition, name)
%ASSERT_TRUE Asserts that condition is true.
    if ~condition
        error('ASSERT_TRUE failed: %s', name);
    end
end