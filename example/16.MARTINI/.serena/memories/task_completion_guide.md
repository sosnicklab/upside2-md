# Task Completion Guide

## When Finishing a Task

### 1. Verify Execution
- Ensure all stages of the simulation completed successfully
- Check for error messages in log files
- Verify output files are created

### 2. Update Documentation
- Update `task_plan.md` to mark completed phases
- Add execution details to `progress.md`
- Add any new findings to `findings.md`

### 3. Code Quality
- Ensure shell scripts have appropriate error handling
- Verify Python scripts follow PEP8 style guidelines
- Check for hardcoded values that should be parameterized

### 4. Test the Workflow
- Run test simulations to verify changes
- Check that input files are correctly generated
- Ensure simulation parameters are properly set

### 5. Performance Check
- Monitor simulation speed and stability
- Check for any unexpected behavior
- Verify that physical properties are reasonable

## Common Post-Task Checks
- Did the simulation reach the target temperature and pressure?
- Are box dimensions changing appropriately?
- Is the system equilibrated before production?
- Are trajectory files properly formatted?
