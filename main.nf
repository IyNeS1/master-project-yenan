process sayHello {
    publishDir 'results', mode: 'copy'
    output:
        path 'output.txt'
    script:
        """
        echo 'Hello World!' > output.txt
        """
}

workflow {
    // Run the sayHello process
    sayHello()
}