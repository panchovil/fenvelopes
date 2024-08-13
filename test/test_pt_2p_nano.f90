program main
    use dtypes, only: envelope
    use constants, only: pr
    use legacy_ar_models, only: nc, z
    use flap, only: command_line_interface
    use stdlib_ansi, only: blue => fg_color_blue, red => fg_color_red, &
                           operator(//), operator(+), &
                           style_reset, style_blink_fast, style_bold, style_underline    
    implicit none
    real(pr) :: et, st

    type(command_line_interface) :: cli
    integer :: cli_error
    character(len=99) :: infile="input.nml"
    character(len=99) :: ouput_path = "test/fenvelopes_nano_output/"
    !character(len=99) ::  ouput_path
    !cli_string
 
    type(envelope) :: pt_bub, pt_dew, pt_hpl !! Shared 2ph-PT envelopes

    real(pr) :: pt_bub_t0 = 180
    real(pr) :: pt_dew_t0 = 180
 
    ! Setup everything
    
    call setup_nano!(infile, ouput_path)
 
    ! PT Envelopes
    call cpu_time(st)
    !call pt_envelopes
    print*, "pt_envelopes"
    call cpu_time(et)
    print *, "PT: ", (et - st)*1000, "cpu ms"
       
    call exit

contains    
subroutine setup_nano!(infile, ouput_path)
    !! Setup system
    !!
    !! Make output folder (if necessary) and/or clean everyhing in an
    !! existing one. Then read input files to setup needed parameters.
    use io_nml, only: read_system, write_system
    !use inj_envelopes, only: setup_inj => from_nml
    integer :: funit_system
    !character(len=99) :: infile, ouput_path

    call system("mkdir -p "//trim(ouput_path))
    call system("rm "//trim(ouput_path)//"*")

    !call setup_cli
    !call cli%get(val=infile, switch="--infile", error=cli_error)

    call read_system(trim(infile))
    !call setup_inj(trim(infile))

    open (newunit=funit_system, file="systemdata.nml")
    call write_system(funit_system)
    close (funit_system)
end subroutine setup_nano
subroutine pt_envelopes_nano
    use envelopes, only: Fnano
end subroutine pt_envelopes_nano
end program main