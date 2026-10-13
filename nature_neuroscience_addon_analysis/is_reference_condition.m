function tf = is_reference_condition(label)
% Reference (expected lower) conditions: Low, SO Peak, Bilateral Trough.
% Used for one-sided directions (e.g. POST < PRE for these conditions).
tf = contains(label, ["Low", "Peak", "Bilateral"]);
end
