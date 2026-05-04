use std::fs;
use std::process::Command;
use std::process::ExitStatus;

fn run_test_case(cmd_args: &str, stdout_file: &str) -> ExitStatus {
    let arg_vec: Vec<String> = cmd_args
        .split_whitespace()
        .map(|s| s.to_string())
        .collect();

    let binary_path = env!("CARGO_BIN_EXE_discoord");

    let output = Command::new(binary_path)
        .args(&arg_vec)
        .output()
        .expect("Failed to run command-line tool");

    fs::write(stdout_file, &output.stdout).expect("Failed to write stdout to file");

    output.status
}

/// Strip all `##`-prefixed header lines so that tests don't care about version
/// strings, thread counts, or other run-time metadata.
fn normalize_log(output: &str) -> String {
    output
        .lines()
        .filter(|line| !line.starts_with("##"))
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn test_single_assembly_fasta() {
    let status = run_test_case(
        "--output-dir tests/data/input/test1/test --map-sequences \
         --log-level detailed \
         --reference-default tests/data/input/ref_dir1/GCF999999.1.fa \
         tests/data/input/test1/input.fa",
        "tests/data/input/test1/stdout_fasta.log",
    );
    assert!(!status.success(), "Command unexpectedly returned success: {:?}", status);

    let expected = fs::read_to_string("tests/data/input/test1/expected-input.fa")
        .expect("Failed to read expected output file");
    let actual = fs::read_to_string("tests/data/input/test1/test/input.fa")
        .expect("Failed to read actual output file");
    assert_eq!(expected, actual, "Output does not match!");

    let expected_log = fs::read_to_string("tests/data/input/test1/expected-detailed.log")
        .expect("Failed to read expected log file");
    let actual_log = fs::read_to_string("tests/data/input/test1/stdout_fasta.log")
        .expect("Failed to read actual log file");
    assert_eq!(normalize_log(&expected_log), normalize_log(&actual_log), "Log does not match!");
}

#[test]
fn test_single_assembly_twobit() {
    let status = run_test_case(
        "--output-dir tests/data/input/test1/test --map-sequences \
         --log-level detailed \
         --reference-default tests/data/input/ref_dir1/GCF999999.1.2bit \
         tests/data/input/test1/input.fa",
        "tests/data/input/test1/stdout_twobit.log",
    );
    assert!(!status.success(), "Command unexpectedly returned success: {:?}", status);

    let expected = fs::read_to_string("tests/data/input/test1/expected-input.fa")
        .expect("Failed to read expected output file");
    let actual = fs::read_to_string("tests/data/input/test1/test/input.fa")
        .expect("Failed to read actual output file");
    assert_eq!(expected, actual, "Output does not match!");

    let expected_log = fs::read_to_string("tests/data/input/test1/expected-detailed.log")
        .expect("Failed to read expected log file");
    let actual_log = fs::read_to_string("tests/data/input/test1/stdout_twobit.log")
        .expect("Failed to read actual log file");
    assert_eq!(normalize_log(&expected_log), normalize_log(&actual_log), "Log does not match!");
}

#[test]
fn test_single_assembly_compressed_fasta() {
    let status = run_test_case(
        "--output-dir tests/data/input/test1/test --map-sequences \
         --log-level detailed \
         --reference-default tests/data/input/ref_dir1/GCF999999.1.2bit \
         tests/data/input/test1/input.fa.gz",
        "tests/data/input/test1/stdout_compressed.log",
    );
    assert!(!status.success(), "Command unexpectedly returned success: {:?}", status);

    let expected = fs::read_to_string("tests/data/input/test1/expected-input.fa")
        .expect("Failed to read expected output file");
    let actual = fs::read_to_string("tests/data/input/test1/test/input.fa")
        .expect("Failed to read actual output file");
    assert_eq!(expected, actual, "Output does not match!");

    let expected_log = fs::read_to_string("tests/data/input/test1/expected-detailed.log")
        .expect("Failed to read expected log file");
    let actual_log = fs::read_to_string("tests/data/input/test1/stdout_compressed.log")
        .expect("Failed to read actual log file");
    assert_eq!(normalize_log(&expected_log), normalize_log(&actual_log), "Log does not match!");
}

#[test]
fn test_single_assembly_stockholm() {
    let status = run_test_case(
        "--output-dir tests/data/input/test2/test --map-sequences \
         --log-level detailed \
         --reference-default tests/data/input/ref_dir1/GCF999999.1.2bit \
         tests/data/input/test2/input.stk",
        "tests/data/input/test2/stdout.log",
    );
    assert!(!status.success(), "Command unexpectedly returned success: {:?}", status);

    let expected = fs::read_to_string("tests/data/input/test2/expected-input.stk")
        .expect("Failed to read expected output file");
    let actual = fs::read_to_string("tests/data/input/test2/test/input.stk")
        .expect("Failed to read actual output file");
    assert_eq!(expected, actual, "Output does not match!");

    let expected_log = fs::read_to_string("tests/data/input/test2/expected-detailed.log")
        .expect("Failed to read expected log file");
    let actual_log = fs::read_to_string("tests/data/input/test2/stdout.log")
        .expect("Failed to read actual log file");
    assert_eq!(normalize_log(&expected_log), normalize_log(&actual_log), "Log does not match!");
}
