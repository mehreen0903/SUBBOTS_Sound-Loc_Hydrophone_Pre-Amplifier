// =============================================================================
// ping_doa.sv  --  High-level framework: 30 kHz underwater ping -> bearing angle
// Target : Terasic DE10-Lite (MAX10 10M50) + TI ADS8528 + 5 hydrophones
//
// ANALOG FRONT END (outside the FPGA, assumed):
//   hydrophone -> preamp -> band-pass (~20-45 kHz) / anti-alias LPF (<~100 kHz)
//   -> ADS8528 input (A0,B0,C0,D0,A1 = 5 hydrophones; unused inputs grounded)
//
// DIGITAL PIPELINE:
//   ads8528_par_if   : one shared CONVST -> all channels sampled SIMULTANEOUSLY,
//                      then read results over the parallel bus (500 kSPS/ch)
//   fir_serial x5    : 64-tap band-pass FIR centered on 30 kHz (1 MAC/channel)
//   ping_detector    : energy/envelope threshold + ring-buffer capture
//   tdoa_xcorr       : cross-correlate ch0 against ch1..4, parabolic sub-sample
//                      interpolation (replaces explicit upsampling)
//   doa_solve        : least-squares plane-wave fit -> (ux,uy) -> CORDIC atan2
//
// NOTE: framework-level code. Verify every ADS8528 timing/pin detail against the
// datasheet, and simulate each block before trusting it.
// =============================================================================

package ping_pkg;
  localparam int NCH      = 5;              // hydrophones
  localparam int CLK_HZ   = 100_000_000;    // from MAX10 PLL (50 MHz board clk x2)
  localparam int FS_HZ    = 500_000;        // per-channel sample rate
  localparam int FIR_TAPS = 64;
  localparam int DEPTH    = 512;            // capture ring-buffer depth / channel
  localparam int AW       = 9;              // log2(DEPTH)
  localparam int XC_W     = 256;            // correlation window (samples)
  localparam int XC_L     = 48;             // max lag +/- (samples) = max baseline / c * Fs
  localparam int PRE_TRIG = 64;             // samples kept before trigger
  localparam int POST     = DEPTH - PRE_TRIG;
  typedef logic signed [15:0] s16_t;
  typedef logic signed [19:0] tau_t;        // TDOA, Q12.8 samples (1 LSB = 1/256 sample)
endpackage


// -----------------------------------------------------------------------------
// ADS8528 parallel interface.
// CONVST A/B/C/D are tied together on the PCB (or driven by this one signal) so
// every channel is sampled at the same instant -> preserves inter-channel phase.
// -----------------------------------------------------------------------------
module ads8528_par_if import ping_pkg::*; (
  input  logic        clk, rst,
  output logic        convst,
  input  logic        busy,
  output logic        cs_n, rd_n,
  input  logic [15:0] db,               // 12-bit result left-justified in DB[15:4]
  output s16_t        sample [NCH],
  output logic        sample_valid
);
  localparam int TICKS = CLK_HZ / FS_HZ;           // 200 clocks @100MHz
  logic [$clog2(TICKS)-1:0] tick_cnt;
  logic tick;
  assign tick = (tick_cnt == TICKS-1);
  always_ff @(posedge clk) tick_cnt <= (rst || tick) ? '0 : tick_cnt + 1'b1;

  logic busy_s1, busy_s2;                          // synchronizer
  always_ff @(posedge clk) {busy_s2, busy_s1} <= {busy_s1, busy};

  typedef enum logic [2:0] {S_IDLE, S_CONV, S_BUSY_HI, S_BUSY_LO, S_RD_LO, S_RD_HI} st_t;
  st_t st;
  logic [2:0] cnt;
  logic [$clog2(NCH)-1:0] idx;

  always_ff @(posedge clk) begin
    sample_valid <= 1'b0;
    if (rst) begin
      st <= S_IDLE; convst <= 1'b0; cs_n <= 1'b1; rd_n <= 1'b1; cnt <= '0; idx <= '0;
    end else case (st)
      S_IDLE:    if (tick) begin convst <= 1'b1; cnt <= '0; st <= S_CONV; end
      S_CONV:    begin                              // hold CONVST high a few clocks
                   cnt <= cnt + 1'b1;
                   if (cnt == 3'd3) begin convst <= 1'b0; st <= S_BUSY_HI; end
                 end
      S_BUSY_HI: if (busy_s2) st <= S_BUSY_LO;      // TODO: timeout if BUSY never rises
      S_BUSY_LO: if (!busy_s2) begin cs_n <= 1'b0; idx <= '0; cnt <= '0; rd_n <= 1'b0; st <= S_RD_LO; end
      S_RD_LO:   begin                              // RD low ~30 ns, then latch
                   cnt <= cnt + 1'b1;
                   if (cnt == 3'd2) begin
                     sample[idx] <= s16_t'($signed(db[15:4]));   // sign-extend 12 -> 16
                     rd_n <= 1'b1; cnt <= '0; st <= S_RD_HI;
                   end
                 end
      S_RD_HI:   begin                              // RD high ~20 ns
                   cnt <= cnt + 1'b1;
                   if (cnt == 3'd1) begin
                     if (idx == NCH-1) begin
                       cs_n <= 1'b1; sample_valid <= 1'b1; st <= S_IDLE;
                       // unread channels (if any) are discarded - check datasheet
                     end else begin
                       idx <= idx + 1'b1; rd_n <= 1'b0; cnt <= '0; st <= S_RD_LO;
                     end
                   end
                 end
      default:   st <= S_IDLE;
    endcase
  end
endmodule


// -----------------------------------------------------------------------------
// Serial (1 MAC) FIR. 200 clocks per sample are available, so 64 taps per channel
// uses only one DSP multiplier per channel.
// Generate coefficients offline (Python: scipy.signal.firwin(64, [25e3,35e3],
// pass_zero=False, fs=500e3)), scale to Q1.15, write as hex to COEF_FILE.
// All 5 channels MUST use identical coefficients (identical group delay).
// -----------------------------------------------------------------------------
module fir_serial import ping_pkg::*; #(
  parameter COEF_FILE = "bp30k_coefs.hex"   // untyped: Quartus rejects "string" here
)(
  input  logic clk, rst,
  input  logic in_valid,
  input  s16_t in_data,
  output logic out_valid,
  output s16_t out_data
);
  logic signed [15:0] coef [0:FIR_TAPS-1];           // coefficient ROM (Q1.15)
  initial $readmemh(COEF_FILE, coef);

  s16_t dl [FIR_TAPS];                               // circular delay line
  logic [$clog2(FIR_TAPS)-1:0] wptr, k;
  logic signed [39:0] acc;
  logic busy;

  logic signed [39:0] acc_next;
  assign acc_next = acc + (dl[wptr - k] * coef[k]);  // wptr-k wraps naturally

  always_ff @(posedge clk) begin
    out_valid <= 1'b0;
    if (rst) begin busy <= 1'b0; wptr <= '0; end
    else if (in_valid && !busy) begin
      dl[wptr] <= in_data; k <= '0; acc <= '0; busy <= 1'b1;
    end else if (busy) begin
      acc <= acc_next;
      if (k == FIR_TAPS-1) begin
        busy      <= 1'b0;
        out_data  <= s16_t'(acc_next >>> 15);        // Q1.15 coefficients
        out_valid <= 1'b1;
        wptr      <= wptr + 1'b1;
      end else k <= k + 1'b1;
    end
  end
endmodule


// -----------------------------------------------------------------------------
// Ping detector + capture memory.
// Envelope = fast-attack / slow-decay of sum(|x|) over all channels (after FIR).
// Ring buffer always records; on trigger it runs POST more samples then freezes.
// The frozen window starts PRE_TRIG samples before the trigger, so the correlation
// window sees the direct-path ping onset (before reverb/multipath).
// -----------------------------------------------------------------------------
module ping_detector import ping_pkg::*; #(
  parameter int THRESH = 4000,        // tune from noise floor; consider adaptive
  parameter int DECAY  = 6,
  parameter int HOLDOFF = 250_000     // samples (0.5 s @500 kSPS): ignore the rest of a
                                      // 4 ms pulse; pinger repeats every 1 s
)(
  input  logic clk, rst,
  input  logic in_valid,
  input  s16_t in_data [NCH],
  input  logic rearm,
  output logic frozen,
  output logic [AW-1:0] start_ptr,
  // two read ports for the correlator (1-cycle latency)
  input  logic [AW-1:0] ref_addr, oth_addr,
  input  logic [$clog2(NCH)-1:0] oth_ch,
  output s16_t ref_data, oth_data
);
  s16_t mem [NCH][DEPTH];
  logic [AW-1:0] wptr;
  logic [19:0] abs_sum, env;
  logic [$clog2(DEPTH)-1:0] post_cnt;
  logic triggered;
  logic [17:0] lock;                                  // hold-off down-counter

  always_comb begin
    abs_sum = '0;
    for (int i = 0; i < NCH; i++)
      abs_sum = abs_sum + 20'(in_data[i][15] ? -in_data[i] : in_data[i]);
  end

  always_ff @(posedge clk) begin
    if (rst) begin
      wptr <= '0; env <= '0; triggered <= 1'b0; frozen <= 1'b0; post_cnt <= '0; lock <= '0;
    end else begin
      if (in_valid && lock != 0) lock <= lock - 1'b1;   // runs even while frozen
      if (rearm) begin frozen <= 1'b0; triggered <= 1'b0; end
      if (in_valid && !frozen) begin
        for (int i = 0; i < NCH; i++) mem[i][wptr] <= in_data[i];
        wptr <= wptr + 1'b1;
        env  <= (abs_sum > env) ? abs_sum : env - (env >> DECAY);
        if (!triggered && lock == 0 && env > THRESH) begin
          triggered <= 1'b1; post_cnt <= '0; lock <= HOLDOFF;
        end
        if (triggered) begin
          post_cnt <= post_cnt + 1'b1;
          if (post_cnt == POST-1) begin frozen <= 1'b1; start_ptr <= wptr + 1'b1; end
        end
      end
    end
  end

  always_ff @(posedge clk) begin                     // read ports
    ref_data <= mem[0][ref_addr];
    oth_data <= mem[oth_ch][oth_addr];
  end
endmodule


// -----------------------------------------------------------------------------
// Small sequential unsigned divider (restoring). Replaces a huge combinational '/'.
// -----------------------------------------------------------------------------
module div_seq #(parameter int N = 56)(
  input  logic clk, rst, start,
  input  logic [N-1:0] num, den,
  output logic [N-1:0] quo,
  output logic done
);
  logic [N-1:0] rem, q;
  logic [$clog2(N+1)-1:0] cnt;
  logic busy;
  logic [N-1:0] r2;
  assign r2 = {rem[N-2:0], q[N-1]};
  always_ff @(posedge clk) begin
    done <= 1'b0;
    if (rst) busy <= 1'b0;
    else if (start && !busy) begin rem <= '0; q <= num; cnt <= N; busy <= 1'b1; end
    else if (busy) begin
      if (r2 >= den) begin rem <= r2 - den; q <= {q[N-2:0], 1'b1}; end
      else           begin rem <= r2;       q <= {q[N-2:0], 1'b0}; end
      cnt <= cnt - 1'b1;
      if (cnt == 1) begin busy <= 1'b0; done <= 1'b1; end
    end
  end
  assign quo = q;
endmodule


// -----------------------------------------------------------------------------
// TDOA by cross-correlation of channel 0 against channels 1..NCH-1.
//   corr[lag] = sum_n ref[n+L] * oth[n+L+lag],  lag = -L..+L
// Integer peak + parabolic interpolation gives ~1/20-1/50 sample resolution
// (at 500 kSPS that is far finer than 2 us), so no explicit upsampler is needed.
// Cost: 4 ch x 97 lags x 256 samples x 2 clk ~ 200k clk = 2 ms @100 MHz.
// Sign: tau[c-1] > 0 means ch c hears the ping LATER than ch 0.
// -----------------------------------------------------------------------------
module tdoa_xcorr import ping_pkg::*; (
  input  logic clk, rst, start,
  input  logic [AW-1:0] start_ptr,
  output logic [AW-1:0] ref_addr, oth_addr,
  output logic [$clog2(NCH)-1:0] oth_ch,
  input  s16_t ref_data, oth_data,
  output tau_t tau [NCH-1],
  output logic [NCH-2:0] tau_ok,
  output logic done
);
  typedef enum logic [3:0] {X_IDLE, X_FETCH, X_MAC, X_STORE, X_PEAK,
                            X_DIVGO, X_DIVWAIT, X_NEXT, X_DONE} st_t;
  st_t st;

  logic signed [39:0] acc;
  logic signed [39:0] corr [2*XC_L+1];
  logic [$clog2(XC_W)-1:0] n;
  int lag;
  int scan_i, pk_idx;
  logic signed [39:0] pk_val;

  logic [$clog2(NCH):0] ch_r;
  assign oth_ch   = ch_r[$clog2(NCH)-1:0];
  assign ref_addr = start_ptr + AW'(XC_L + n);
  assign oth_addr = start_ptr + AW'(XC_L + n + lag);   // modulo arithmetic handles negative lag

  // interpolation
  logic signed [41:0] cm, c0, cp;
  logic signed [42:0] num_s, den_s;
  logic        div_start, div_done, neg;
  logic [55:0] div_num, div_den, div_quo;
  div_seq #(.N(56)) u_div (.clk(clk), .rst(rst), .start(div_start), .num(div_num),
                           .den(div_den), .quo(div_quo), .done(div_done));

  always_ff @(posedge clk) begin
    div_start <= 1'b0; done <= 1'b0;
    if (rst) st <= X_IDLE;
    else case (st)
      X_IDLE:  if (start) begin
                 ch_r <= 1; lag <= -XC_L; n <= '0; acc <= '0; tau_ok <= '0; st <= X_FETCH;
               end
      X_FETCH: st <= X_MAC;                                // BRAM read latency
      X_MAC:   begin
                 acc <= acc + (ref_data * oth_data);
                 if (n == XC_W-1) st <= X_STORE; else begin n <= n + 1'b1; st <= X_FETCH; end
               end
      X_STORE: begin
                 corr[lag + XC_L] <= acc;
                 if (lag == XC_L) begin
                   pk_val <= 40'sh80_0000_0000; pk_idx <= 0; scan_i <= 0; st <= X_PEAK;
                 end else begin lag <= lag + 1; n <= '0; acc <= '0; st <= X_FETCH; end
               end
      X_PEAK:  begin
                 if (corr[scan_i] > pk_val) begin pk_val <= corr[scan_i]; pk_idx <= scan_i; end
                 if (scan_i == 2*XC_L) st <= X_DIVGO; else scan_i <= scan_i + 1;
               end
      X_DIVGO: begin
                 if (pk_idx == 0 || pk_idx == 2*XC_L) st <= X_NEXT;     // peak at edge: invalid
                 else begin
                   cm = corr[pk_idx-1]; c0 = corr[pk_idx]; cp = corr[pk_idx+1];
                   num_s = cm - cp;                                     // delta = (cm-cp) / (2(cm-2c0+cp))
                   den_s = 2 * (cm - 2*c0 + cp);
                   neg   <= num_s[42] ^ den_s[42];
                   div_num   <= 56'((num_s[42] ? -num_s : num_s)) << 8;  // Q8 fraction
                   div_den   <= 56'(den_s[42] ? -den_s : den_s);
                   if (den_s == 0) st <= X_NEXT;
                   else begin div_start <= 1'b1; st <= X_DIVWAIT; end
                 end
               end
      X_DIVWAIT: if (div_done) begin
                 tau[ch_r-1]    <= tau_t'(((pk_idx - XC_L) <<< 8) + (neg ? -$signed({1'b0, div_quo[19:0]})
                                                                         :  $signed({1'b0, div_quo[19:0]})));
                 tau_ok[ch_r-1] <= 1'b1;
                 st <= X_NEXT;
               end
      X_NEXT:  begin
                 if (ch_r == NCH-1) st <= X_DONE;
                 else begin ch_r <= ch_r + 1'b1; lag <= -XC_L; n <= '0; acc <= '0; st <= X_FETCH; end
               end
      X_DONE:  begin done <= 1'b1; st <= X_IDLE; end
      default: st <= X_IDLE;
    endcase
  end
endmodule


// -----------------------------------------------------------------------------
// Iterative CORDIC atan2 (vectoring mode). Output: 16-bit full-circle angle,
// 0..65535 = 0..360 deg. Only direction matters, so input scaling is arbitrary.
// -----------------------------------------------------------------------------
module cordic_atan2 #(parameter int W = 42, ITER = 16)(
  input  logic clk, rst, start,
  input  logic signed [W-1:0] x_in, y_in,
  output logic [15:0] angle,
  output logic done
);
  // round(atan(2^-i) / 2pi * 65536)
  localparam logic [15:0] ATAN [16] = '{16'd8192,16'd4836,16'd2555,16'd1297,16'd651,16'd326,
                                        16'd163,16'd81,16'd41,16'd20,16'd10,16'd5,16'd3,16'd1,16'd1,16'd0};
  logic signed [W-1:0] x, y;
  logic [15:0] z;
  logic [4:0] i;
  logic busy;
  always_ff @(posedge clk) begin
    done <= 1'b0;
    if (rst) busy <= 1'b0;
    else if (start && !busy) begin
      // pre-rotate by 180 deg if x < 0 so iterations converge
      if (x_in < 0) begin x <= -x_in; y <= -y_in; z <= 16'd32768; end
      else          begin x <=  x_in; y <=  y_in; z <= 16'd0;     end
      i <= '0; busy <= 1'b1;
    end else if (busy) begin
      if (y >= 0) begin x <= x + (y >>> i); y <= y - (x >>> i); z <= z + ATAN[i]; end
      else        begin x <= x - (y >>> i); y <= y + (x >>> i); z <= z - ATAN[i]; end
      i <= i + 1'b1;
      if (i == ITER-1) begin busy <= 1'b0; angle <= z; done <= 1'b1; end   // 1 clk lag on last z: negligible
    end
  end
endmodule


// -----------------------------------------------------------------------------
// Direction solve. Plane wave from unit vector u (pointing TO the source):
//   tau_i = -( (p_i - p_0) . u ) / c
// => u = G * tau, with G = -c*Ts * pinv(D)   (2 x 4, D rows = p_i - p_0)
// Generate G offline from the real hydrophone coordinates (Python/numpy), Q2.14.
// For 3-D bearing (azimuth+elevation) make G 3x4 and add a second atan2 on
// (uz, hypot(ux,uy)).
// -----------------------------------------------------------------------------
module doa_solve import ping_pkg::*; (
  input  logic clk, rst,
  input  logic in_valid,
  input  tau_t tau [NCH-1],
  output logic [15:0] angle,        // 0..65535 = 0..360 deg, 0 = array +x axis
  output logic angle_valid
);
  // *** PLACEHOLDER - REPLACE with values computed from your array geometry ***
  localparam logic signed [15:0] G [2][NCH-1] = '{
    '{ 16'sd4096,  16'sd0,    -16'sd4096,  16'sd0    },
    '{ 16'sd0,     16'sd4096,  16'sd0,    -16'sd4096 }
  };

  logic signed [41:0] ux, uy;
  logic signed [41:0] ux_r, uy_r;
  logic cordic_start;
  always_ff @(posedge clk) begin
    cordic_start <= 1'b0;
    if (in_valid) begin
      ux = '0; uy = '0;
      for (int i = 0; i < NCH-1; i++) begin
        ux = ux + G[0][i] * tau[i];
        uy = uy + G[1][i] * tau[i];
      end
      cordic_start <= 1'b1;
      ux_r <= ux; uy_r <= uy;
    end
  end
  cordic_atan2 #(.W(42)) u_cordic (.clk(clk), .rst(rst), .start(cordic_start),
                                   .x_in(ux_r), .y_in(uy_r), .angle(angle), .done(angle_valid));
endmodule


// -----------------------------------------------------------------------------
// Top level (DE10-Lite). Map ports to GPIO header pins in the .qsf.
// ADS8528 static pins (tie on PCB/jumpers): parallel mode select, hardware mode,
// internal reference enable, input range, CONVST A-D tied together.
// -----------------------------------------------------------------------------
module ping_doa_top import ping_pkg::*; (
  input  logic        MAX10_CLK1_50,
  input  logic [1:0]  KEY,                // KEY[0] = reset (active low)
  output logic [9:0]  LEDR,
  // ADS8528
  output logic        ADC_CONVST,
  input  logic        ADC_BUSY,
  output logic        ADC_CS_N, ADC_RD_N,
  output logic        ADC_WR_N,           // tie high in hardware mode
  output logic        ADC_RESET,
  input  logic [15:0] ADC_DB,
  // result (also route to UART/GPIO as needed)
  output logic [15:0] ANGLE_OUT,
  output logic        ANGLE_STROBE
);
  // 100 MHz from the MAX10 PLL (generate "pll_100m" with the IP Catalog)
  logic clk, pll_locked;
  pll_100m u_pll (.inclk0(MAX10_CLK1_50), .c0(clk), .locked(pll_locked));

  logic [3:0] rst_sr;
  logic rst;
  always_ff @(posedge clk or negedge KEY[0])
    if (!KEY[0]) rst_sr <= '1; else rst_sr <= {rst_sr[2:0], ~pll_locked};
  assign rst       = rst_sr[3];
  assign ADC_RESET = rst;                 // ADS8528 RESET is active high
  assign ADC_WR_N  = 1'b1;

  // --- acquisition ---
  s16_t raw [NCH];  logic raw_valid;
  ads8528_par_if u_adc (.clk(clk), .rst(rst), .convst(ADC_CONVST), .busy(ADC_BUSY),
                        .cs_n(ADC_CS_N), .rd_n(ADC_RD_N), .db(ADC_DB),
                        .sample(raw), .sample_valid(raw_valid));

  // --- band-pass FIR per channel ---
  s16_t filt [NCH];  logic [NCH-1:0] filt_v;
  genvar c;
  generate
    for (c = 0; c < NCH; c = c + 1) begin : g_fir
      fir_serial #(.COEF_FILE("bp30k_coefs.hex")) u_fir (
        .clk(clk), .rst(rst), .in_valid(raw_valid), .in_data(raw[c]),
        .out_valid(filt_v[c]), .out_data(filt[c]));
    end
  endgenerate

  // --- detect + capture ---
  logic frozen, rearm;
  logic [AW-1:0] start_ptr, ref_addr, oth_addr;
  logic [$clog2(NCH)-1:0] oth_ch;
  s16_t ref_data, oth_data;
  ping_detector u_det (.clk(clk), .rst(rst), .in_valid(filt_v[0]), .in_data(filt), .rearm(rearm),
                       .frozen(frozen), .start_ptr(start_ptr), .ref_addr(ref_addr),
                       .oth_addr(oth_addr), .oth_ch(oth_ch),
                       .ref_data(ref_data), .oth_data(oth_data));

  // --- TDOA ---
  logic xc_start, xc_done, frozen_d;
  tau_t tau [NCH-1];  logic [NCH-2:0] tau_ok;
  always_ff @(posedge clk) begin frozen_d <= frozen; xc_start <= frozen & ~frozen_d; end
  tdoa_xcorr u_xc (.clk(clk), .rst(rst), .start(xc_start), .start_ptr(start_ptr),
                   .ref_addr(ref_addr), .oth_addr(oth_addr), .oth_ch(oth_ch),
                   .ref_data(ref_data), .oth_data(oth_data), .tau(tau),
                   .tau_ok(tau_ok), .done(xc_done));

  // --- DOA (only if all TDOAs valid) ---
  logic doa_go;
  assign doa_go = xc_done & (&tau_ok);
  doa_solve u_doa (.clk(clk), .rst(rst), .in_valid(doa_go), .tau(tau), .angle(ANGLE_OUT), .angle_valid(ANGLE_STROBE));

  // re-arm detector once the solver has consumed the capture
  assign rearm = xc_done;

  // simple debug: show top angle bits on LEDs
  always_ff @(posedge clk) if (ANGLE_STROBE) LEDR <= ANGLE_OUT[15:6];
endmodule